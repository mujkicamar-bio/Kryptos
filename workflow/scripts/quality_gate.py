"""S5: the quality gate. Halts the run if the cascade has a recall problem.

POSITIVE CONTROL (success criterion SC2)
    Known plasmid biology is run through the same cascade as everything else and must come
    out FUNCTIONAL. ECLIPSE recovered 99.2-100% of 246 virulence, 42 AMR and 75 essential
    genes this way. Anything KNOWN that emerges dark is a recall failure, and a target list
    built on a broken annotation step is worse than no target list - it costs bench time
    and it is not detectably wrong until the assays fail.

    Its absence from v1 was the single point of unanimous reviewer criticism.

WHY THERE IS NO LONGER A BACKBONE STOP-LIST

An earlier version also excluded any protein carrying a curated plasmid-backbone family
name, as a second guard behind the classification. It has been removed deliberately.

The pipeline is two steps: annotate everything that can be annotated, then report against
criteria. If step one works, a complete replication initiator is FUNCTIONAL and is simply
not dark - it needs no list to keep it out. If step one is broken, a hand-maintained list of
family names is the wrong repair, and a dangerous one: an exclusion by name removes the
protein from the record entirely, so its removal is invisible and unauditable. That is not
theoretical - nine of the seventy-three names on that list did not exist in Pfam-A at all,
so those rules had never fired, and nobody could tell from any output.

The curated family table is gone. Pfam-A 38.2 holds 30,134 families, of which 67 mention
replication in their description and 42 mention conjugation; the list named 16 and 15, and
nine of its 73 names did not exist in Pfam-A at all. Functional labels now come from the
tools themselves, in results/s4c/protein_labels.tsv, and are grouped into categories from
the observed vocabulary rather than by hand - and they LABEL a protein that stays in the
table rather than excluding it. Nothing excludes anything by name.

Nothing is deleted here (P5). Proteins are flagged, and the flags are counted.
"""
import _ctx  # noqa: F401
import collections
import csv

from plasmidann.controls import control_recall

cfg = snakemake.params.gate

# ------------------------------------------------------------------------------------
# Read the cascade output.
# ------------------------------------------------------------------------------------
rows = []
with open(snakemake.input.prot, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rows.append(r)

# ------------------------------------------------------------------------------------
# Positive control: reviewed Swiss-Prot proteins of known function, spiked into the query
# set at S2c and carried through every tier exactly as a real protein is.
#
# Independent by construction. An earlier version used proteins the cascade itself had
# labelled as backbone, which could not detect the failure that matters - a protein the
# cascade MISSED never enters a self-drawn control set.
# ------------------------------------------------------------------------------------
CONTROL_PREFIX = "CTRL_"
control_rows = [r for r in rows if r["seq_id"].startswith(CONTROL_PREFIX)]

if not control_rows and cfg.get("require_control_set", True):
    raise SystemExit(
        "S5 QUALITY GATE: no control proteins found in the annotation table.\n"
        "The controls are spiked in at S2c and must reach S3. A gate that silently skips "
        "itself is not a gate - check that results/s2/cascade_input.faa was used as the "
        "cascade query set.")

recall = control_recall(control_rows) if control_rows else 0.0
failed = [r["seq_id"] for r in control_rows if r["functional_class"] != "FUNCTIONAL"]

# WHICH TIER resolved each control, not just whether one did.
#
# The controls are Swiss-Prot proteins and T3 searches swissprot.dmnd, so they self-hit
# there at essentially perfect identity. That is a real annotation path and not a bug, but
# it makes overall recall a weak test: the gate would pass even if the curated Pfam tiers
# were completely broken, because T3 would rescue every control on its own.
#
# Reporting resolution per tier turns a nearly-trivial pass into a diagnostic. If controls
# only ever resolve at T3, T1 and T2 have a recall problem that overall recall hides.
by_tier = collections.Counter(r.get("annot_tier") or "UNRESOLVED" for r in control_rows)
shallow = sum(n for t, n in by_tier.items() if t in ("T1", "T2"))
shallow_fraction = round(shallow / len(control_rows), 4) if control_rows else 0.0

# ------------------------------------------------------------------------------------
# Artefact flags from S2b and edge-partial ORFs are the other two exclusions.
# ------------------------------------------------------------------------------------
artefact_ids = set()
with open(snakemake.input.artefact, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["artefact_flag"] == "1":
            artefact_ids.add(r["seq_id"])

cols = ["seq_id", "is_artefact", "target_eligible", "exclusion_reason"]
n_eligible = 0
with open(snakemake.output.flags, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in rows:
        sid = r["seq_id"]
        if sid.startswith(CONTROL_PREFIX):
            continue          # instrumentation, not a candidate
        reasons = []
        if sid in artefact_ids:
            reasons.append("artefact")
        # Only proteins nothing could name are screening candidates at all.
        if r["functional_class"] not in ("UNCHARACTERIZED_HOMOLOG", "NONE"):
            reasons.append("annotated")
        eligible = int(not reasons)
        n_eligible += eligible
        w.writerow({"seq_id": sid,
                    "is_artefact": int(sid in artefact_ids),
                    "target_eligible": eligible,
                    "exclusion_reason": ",".join(reasons)})

summary = (f"proteins={len(rows)} artefact={len(artefact_ids)} "
           f"target_eligible={n_eligible} control_n={len(control_rows)} "
           f"control_recall={recall}")
with open(snakemake.output.report, "w") as out:
    out.write(summary + "\n")
    out.write(f"min_control_recall={cfg['min_control_recall']}\n")
    out.write(f"\ncontrols resolved per tier (T3 is swissprot, where controls self-hit):\n")
    for tier, n in sorted(by_tier.items()):
        out.write(f"  {tier}\t{n}\n")
    out.write(f"resolved at T1 or T2 (curated Pfam path)\t{shallow}/{len(control_rows)}"
              f"\t{shallow_fraction}\n")
    if failed:
        out.write(f"\ncontrols NOT classed FUNCTIONAL ({len(failed)}):\n")
        for sid in failed[:50]:
            out.write(f"  {sid}\n")
print(summary)

# The gate. Failing here stops the run rather than producing a target list nobody should
# act on.
if recall < cfg["min_control_recall"]:
    raise SystemExit(
        f"S5 QUALITY GATE FAILED: positive control recall {recall} is below the required "
        f"{cfg['min_control_recall']}.\n\n"
        f"{len(failed)} of {len(control_rows)} reviewed Swiss-Prot proteins of KNOWN "
        "function were not classed FUNCTIONAL by this cascade. That means annotation is "
        "losing proteins it should catch, so 'unannotated' carries no information and any "
        "target list built on it would send unannotatable noise to the bench.\n\n"
        f"First failures: {', '.join(failed[:10])}\n\n"
        "NOTE: this gate is about proteins that SHOULD have been annotated. Genuinely "
        "novel proteins - integron cassettes, defence accessories, small ORFs - are "
        "supposed to come out unannotated and are never a reason to halt.\n"
        f"See {snakemake.output.report} for the full list.")
