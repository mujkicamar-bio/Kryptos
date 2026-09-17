"""S9: assemble the evidence for every dark family and select the screening portfolio.

NO SCORES. NO WEIGHTS.

An earlier version reduced four evidence groups to a weighted composite. Every coefficient
in it was invented - 0.40 for selection, 0.25 for lineage breadth, caps at 3 MOB clusters,
10 members, 400 residues - and none was measured, defensible or falsifiable. It also
destroyed the information a screening decision needs: two candidates scoring 0.6 can be
completely different bets, one with overwhelming evidence of being real and no idea what it
does, the other with a sharp hypothesis resting on almost nothing. Those want different
experiments.

What each family carries instead:

  reality_n / reality_lines   how many of four independent, named tests fired, and which
  darkness_state              DARK_NO_STRUCTURE or DARK_FOLD_KNOWN
  hypothesis + conservation   the actual hypothesis and how well it holds across members
  liabilities                 a list of expression risks, not a number
  stratum                     what the microscope is expected to see

Selection is eligibility (a declared count of evidence lines), then a lexicographic ranking
within each stratum. Sensitivity is reported on the one declared integer that remains.
"""
import _ctx  # noqa: F401
import collections
import csv

from plasmidann.targets import (reality_lines, darkness_state, hypothesis_confidence,
                                select_portfolio, eligibility_sensitivity, RANK_PRIORITY,
                                REALITY_TESTS, binds_nucleic_acid, reality_thresholds,
                                check_reality_config)
from plasmidann.peptide import (net_charge, gravy, max_hydrophobic_window, has_tm_helix,
                                is_cationic_amphipathic)

cfg = snakemake.params.prioritisation
port = snakemake.params.portfolio
ctx_cfg = snakemake.params.context
lib_cfg = snakemake.params.library
evo_cfg = snakemake.params.evolution

# Every reality-test threshold, assembled from the config blocks that own them. Two of the
# three come from the evolution block rather than being restated, so the nesting the
# eligibility count relies on cannot drift apart from what S7b actually measured.
check_reality_config(cfg, evo_cfg)
REALITY_THRESHOLDS = reality_thresholds(cfg, evo_cfg)


def read_rows(path):
    """Every row of a TSV, in file order."""
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_tsv(path, key):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r[key]] = r
    return out


families = read_tsv(snakemake.input.families, "family_id")
evolution = read_tsv(snakemake.input.evolution, "family_id")
def read_context(path, ctx_cfg):
    """Reduce the LONG family_context table to one row per family.

    S8c writes one row per (family, category), because the label vocabulary is open and
    cannot be columns. The strongest association per family is selected HERE rather than
    there, because "strongest" depends on the q-value threshold, which is a prioritisation
    parameter and is swept.

    Selection is by q-value, not by conservation or enrichment. A category conserved across
    three plasmids and one conserved across three hundred both read 1.0 conserved, and the
    q-value is the only one of the three numbers that tells them apart.
    """
    rows = collections.defaultdict(dict)
    for r in read_rows(path):
        family = r["family_id"]
        category = r["category"]
        observed = float(r["observed_rate"] or 0.0)
        # Every category the family carries is kept as cons_<category>, so the existing
        # cassette tests still read cons_integron and the report still reads cons_defence.
        rows[family][f"cons_{category}"] = observed
        rows[family][f"enrich_{category}"] = r["enrichment"]
        rows[family][f"q_{category}"] = r["q_value"]
        rows[family][f"sub_{category}"] = r["subcategories"]

        if r["status"] != "SUCCESS":
            continue
        try:
            q = float(r["q_value"])
        except (TypeError, ValueError):
            continue
        if q > ctx_cfg["max_q_value"]:
            continue
        if observed < ctx_cfg["min_context_conservation"]:
            continue
        best = rows[family].get("_best_q")
        if best is None or q < best:
            rows[family]["_best_q"] = q
            rows[family]["top_hypothesis"] = category
            rows[family]["top_conservation"] = observed
            rows[family]["top_enrichment"] = r["enrichment"]
            rows[family]["top_q_value"] = r["q_value"]
            rows[family]["top_subcategories"] = r["subcategories"]
    return rows


context = read_context(snakemake.input.context, ctx_cfg)
structure = read_tsv(snakemake.input.structure, "seq_id")

sequences, name, buf = {}, None, []
for line in open(snakemake.input.faa):
    if line[0] == ">":
        if name:
            sequences[name] = "".join(buf)
        name, buf = line[1:].split()[0], []
    else:
        buf.append(line.strip())
if name:
    sequences[name] = "".join(buf)


def _f(d, k, default=0.0):
    try:
        return float(d.get(k, "") or default)
    except (ValueError, TypeError):
        return default


candidates = []
for fid, fam in families.items():
    rep = fam["representative"]
    seq = sequences.get(rep, "")
    if not seq:
        continue
    evo = evolution.get(fid, {})
    ctx = context.get(fid, {})
    struct = structure.get(rep, {})

    # ---- the record the evidence tests read ----------------------------------------
    record = {
        "dnds_status": evo.get("dnds_status", ""),
        "dnds_median": evo.get("dnds_median", ""),
        "n_mob_clusters": fam["n_mob_clusters"],
        "n_members": fam["n_members"],
        "structural_match": struct.get("target", ""),
    }
    n_reality, fired, implied = reality_lines(record, REALITY_THRESHOLDS)

    # ---- hypothesis: the thing itself, not a number ---------------------------------
    top = ctx.get("top_hypothesis", "")
    top_cons = _f(ctx, "top_conservation")

    # TWO DIFFERENT CLAIMS ABOUT CASSETTES, deliberately kept apart.
    #
    #   in_cassette        at least one member sits in a cassette array. A LABEL, kept
    #                      because a reader explicitly wants to see that a dark ORF is in a
    #                      cassette, and never a filter.
    #   cassette_conserved the association holds across the family, at the same
    #                      conservation bar every other context claim in this pipeline has
    #                      to clear. This is what may decide a stratum.
    #
    # They were one test, `cons_integron > 0`, checked second in the stratum chain - so a
    # single coincidental member out of ten diverted the whole family into the
    # integron_cassette stratum, 175 of the 1,000 constructs.
    cons_integron = _f(ctx, "cons_integron")
    in_cassette = cons_integron > 0
    cassette_conserved = cons_integron >= ctx_cfg["min_context_conservation"]

    if not top and struct.get("target"):
        # A recognisable fold is a hypothesis even with no genomic context - but only if it
        # says what the fold IS. `structural:12as-assembly1_A` is a PDB accession and tells
        # a bench scientist nothing; the description is the hypothesis.
        top = f"structural:{struct.get('target_description') or struct['target']}"
        top_cons = 0.0
    confidence = hypothesis_confidence(top_cons, ctx_cfg) if ctx.get("top_hypothesis") \
        else ("MODERATE" if top else "NONE")

    # ---- liabilities: flags a wet-lab reader acts on, not a scalar ------------------
    length = len(seq)
    tm = has_tm_helix(seq)
    amp = is_cationic_amphipathic(seq)
    liabilities = []
    if length > lib_cfg["length_liability_above_aa"]:
        liabilities.append("long")
    if tm:
        liabilities.append("membrane_may_need_detergent")
    if seq.count("C") >= 6:
        liabilities.append("cysteine_rich")
    if not top:
        liabilities.append("no_hypothesis")

    # ---- stratum: what the microscope is expected to see ----------------------------
    # A DECLARED PRIORITY ORDER, and every match is recorded, not just the winner.
    #
    # A candidate can satisfy several of these at once, and a first-match chain that threw
    # the rest away had a specific failure mode: has_tm_helix is an explicitly coarse
    # hydrophobicity heuristic that over-fires on soluble proteins with a hydrophobic core,
    # and it is checked third - so a false TM call could starve defence_island,
    # nucleic_acid_binding and novel_fold with no trace of what the candidate would
    # otherwise have been. The order still decides; the alternatives now survive into the
    # table so the decision is auditable and a stratum can be re-cut without a re-run.
    #
    # Cassette membership sits near the top because it is the strongest plasmid-specific
    # evidence available - an attC site means the gene was excised, mobilised and retained.
    # It is never a reason to exclude; it only decides which plate the protein goes on.
    matched = []
    if amp and length <= 100:
        matched.append("small_cationic_peptide")
    if ctx.get("top_hypothesis") == "integron" or cassette_conserved:
        matched.append("integron_cassette")
    if tm:
        matched.append("membrane_or_secreted")
    # The toxin_or_ta_adjacent stratum was filled by a 'ta_candidate' hypothesis that
    # fired when a dark ORF's two-gene partner appeared on a hand-written list of 73 Pfam
    # family names. That list is deleted, and its replacement - a toxin_antitoxin CATEGORY
    # derived from the observed vocabulary in results/s4c/protein_labels.tsv - does not
    # exist until a full annotation run has been reviewed. The stratum is declared
    # suspended in config/targets.yaml rather than being filled by a rule with no
    # definition behind it. Restoring it needs only the category, not a code change:
    # `if ctx.get("top_hypothesis") == "toxin_antitoxin"`.
    if ctx.get("top_hypothesis") == "defence":
        matched.append("defence_island")
    # From the DESCRIPTION, never the accession. Testing `"nucle" in target` asked whether a
    # PDB accession such as `200l-assembly1_A` contained the word, which none does, so this
    # stratum - 75 constructs - could never be filled by anything.
    if binds_nucleic_acid(struct.get("target_description")):
        matched.append("nucleic_acid_binding")
    if not struct.get("target"):
        matched.append("novel_fold")
    # A recognised fold that is none of the above: no phenotype prediction, but real
    # structural evidence. Carries no quota, is reported, and can fill a shortfall.
    stratum = matched[0] if matched else "unassigned"

    candidates.append({
        "id": fid, "stratum": stratum,
        "stratum_alternatives": ",".join(matched[1:]) or "none",
        "representative": rep,
        # reality_n counts only INDEPENDENT lines, which is what min_reality_lines means.
        # reality_lines_implied records what else fired but added nothing, so the table
        # still shows everything that was true about the family.
        "reality_n": n_reality, "reality_lines": "+".join(fired) or "none",
        "reality_lines_implied": "+".join(implied) or "none",
        "darkness_state": darkness_state(record),
        "hypothesis": top, "hypothesis_conservation": top_cons,
        "hypothesis_enrichment": ctx.get(f"enrich_{ctx.get('top_hypothesis', '')}", ""),
        "hypothesis_confidence": confidence,
        "hyp_rank": {"HIGH": 2, "MODERATE": 1, "NONE": 0}[confidence],
        "liabilities": ",".join(liabilities) or "none", "liability_n": len(liabilities),
        "length": length, "n_members": int(fam["n_members"]),
        # Gene copies behind the cluster, which is not the member count: one unique
        # sequence carried by 200 plasmids is n_members=1, n_orfs=200. Without both, a
        # widely carried conserved protein is indistinguishable from a singleton.
        "n_orfs": int(fam.get("n_orfs") or fam["n_members"]),
        "n_mob_clusters": int(fam["n_mob_clusters"]),
        "dnds_median": evo.get("dnds_median", ""), "dnds_status": evo.get("dnds_status", ""),
        "in_integron_cassette": int(in_cassette),
        "cassette_conserved": int(cassette_conserved),
        "cons_integron": _f(ctx, "cons_integron"), "cons_defence": _f(ctx, "cons_defence"),
        "structural_match": struct.get("target", ""),
        "structural_description": struct.get("target_description", ""),
        "structure_evalue": struct.get("evalue", ""),
        "net_charge": net_charge(seq), "gravy": gravy(seq),
        "max_hydrophobic_window": max_hydrophobic_window(seq),
        "predicted_tm": int(tm), "topology_method": "hydrophobicity_window",
    })

picked, report = select_portfolio(
    candidates, port["strata"],
    min_reality_lines=cfg["min_reality_lines"],
    reallocate=port.get("reallocate_shortfall", False))
picked_ids = {c["id"] for c in picked}

# The honest replacement for sweeping invented weights. One declared integer, and it means
# something a reader can reason about: how many independent lines of evidence are demanded.
sensitivity = eligibility_sensitivity(candidates, port["strata"])

cols = ["family_id", "representative", "selected", "stratum", "stratum_alternatives",
        "reality_n", "reality_lines", "reality_lines_implied", "darkness_state",
        "hypothesis", "hypothesis_confidence", "hypothesis_conservation",
        "hypothesis_enrichment", "liabilities", "liability_n", "is_pareto_optimal",
        "length", "n_members", "n_orfs", "n_mob_clusters", "dnds_median", "dnds_status",
        "in_integron_cassette", "cassette_conserved", "cons_integron", "cons_defence",
        "structural_match", "structural_description", "structure_evalue",
        "net_charge", "gravy",
        "max_hydrophobic_window", "predicted_tm", "topology_method",
        "min_reality_lines_applied"]

with open(snakemake.output.scored, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for c in sorted(candidates, key=lambda x: (-x["reality_n"],
                                               -x["hypothesis_conservation"])):
        c["selected"] = int(c["id"] in picked_ids)
        c["family_id"] = c["id"]
        c.setdefault("is_pareto_optimal", 0)
        c["min_reality_lines_applied"] = cfg["min_reality_lines"]
        w.writerow(c)

counts = collections.Counter(c["stratum"] for c in picked)
with open(snakemake.output.report, "w") as out:
    out.write("SELECTION PROCEDURE\n")
    out.write(f"  eligibility     >= {cfg['min_reality_lines']} of "
              f"{len(REALITY_TESTS)} independent evidence lines\n")
    out.write(f"  ranking         lexicographic: {' > '.join(RANK_PRIORITY)}\n")
    out.write("  weights         none\n\n")
    out.write(f"candidate families\t{len(candidates)}\n")
    out.write(f"selected\t{len(picked)}\ttarget\t{port['total']}\n")
    out.write(f"reallocate_shortfall\t{port.get('reallocate_shortfall', False)}\n\n")
    out.write("stratum\tselected\tquota\tshortfall\tavailable\n")
    for s, r in report.items():
        out.write(f"{s}\t{r['selected']}\t{r['quota']}\t{r['shortfall']}\t{r['available']}\n")
    out.write("\nEVIDENCE LINES ACROSS SELECTED CANDIDATES\n")
    line_counts = collections.Counter()
    for c in picked:
        for name in c["reality_lines"].split("+"):
            line_counts[name] += 1
    for name, n in line_counts.most_common():
        out.write(f"  {name}\t{n}\n")
    out.write("\nSENSITIVITY TO THE ELIGIBILITY FLOOR\n")
    out.write("  lines demanded\tcandidates selectable\n")
    for n, k in sorted(sensitivity.items()):
        out.write(f"  {n}\t{k}\n")
    # A protein whose sequence is identical across many plasmids is ONE cluster member: it
    # is ORPHAN, fails is_family, and cannot have a dN/dS measured because there is nothing
    # to align. It can therefore reach at most one reality line - multi_lineage - and is
    # ineligible at any floor above 1, despite near-universal conservation being among the
    # strongest prevalence signals in the collection.
    #
    # That is arguably the right definition of a family and definitely the wrong silence.
    # The count is reported so a reader can see the size of what the definition costs.
    widely_carried = [c for c in candidates
                      if c["n_members"] < REALITY_THRESHOLDS["min_members"]
                      and c["n_orfs"] >= 10]
    excluded = [c for c in widely_carried if c["id"] not in picked_ids]
    out.write("\nWIDELY CARRIED BUT NOT A FAMILY\n")
    out.write("  A single unique sequence on many plasmids clusters alone, so it fails\n"
              "  is_family and cannot have dN/dS measured. Reported because the definition\n"
              "  excludes exactly the most conserved proteins in the collection.\n")
    out.write(f"  candidates with <{REALITY_THRESHOLDS['min_members']} members "
              f"but >=10 ORF copies\t{len(widely_carried)}\n")
    out.write(f"  of those, not selected\t{len(excluded)}\n")

    n_cass = sum(1 for c in candidates if c["in_integron_cassette"])
    out.write(f"\ndark families in an integron cassette\t{n_cass}\n")
    out.write(f"  of which selected\t{sum(1 for c in picked if c['in_integron_cassette'])}\n")

print(f"families={len(candidates)} selected={len(picked)}/{port['total']} "
      f"(eligibility >= {cfg['min_reality_lines']} evidence lines, no weights)")
for s, r in report.items():
    print(f"  {s:<26} {r['selected']:>5} / {r['quota']:<5} available={r['available']}")
