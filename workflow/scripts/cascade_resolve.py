"""S3 final step: turn every tier's hits into one row per protein.

This is where the cascade's evidence becomes a decision. Two principles govern it:

  * The class is a property of the PROTEIN. Coverage is merged across all informative hits
    from all tiers before anything is decided. Classifying per hit put 17% of Pfam-hit
    proteins into DOMAIN_ONLY at >=50% explained - and since DOMAIN_ONLY is
    target-eligible, complete replication initiators and T4SS ATPases were walking into
    the screening pool.

  * min_explained is applied HERE, post hoc, not during the search. The search narrowed on
    narrow_at (0.7), so every protein explained below 0.7 was seen by every tier and this
    threshold can be swept up to 0.7 without re-running anything.

Two columns describe coverage, over disjoint evidence, and they are not redundant:

  explained_fraction   merged INFORMATIVE spans / length. Informative about ANNOTATED
                       proteins. Constant 0 across the dark set by construction.
  dark_covered_fraction merged UNINFORMATIVE spans / length. The discriminating axis for
                       the dark set: 95% covered by 'hypothetical protein' across three
                       databases is a real, conserved, full-length unnamed protein; a
                       single 20-aa fragment hit is not.

EVERY UNIQUE PROTEIN GETS A ROW, and annot_source says where it came from:

  self            searched by the cascade (a search representative, or a control)
  representative  a member of a 90% search cluster (S2s); the row is its representative's,
                  named in annot_representative. Coverage fields describe the
                  representative, which is within ~20% of the member's length (cov-mode 0)
  plasmidscope    Tier 0
  not_searched    outside every family with an unexplained small-plasmid protein;
                  functional_class NOT_SEARCHED, which is neither dark nor annotated
  artefact_antifam  flagged by AntiFam (S2b), so it skipped every annotation tier, Tier 0
                  included; functional_class NOT_SEARCHED
"""
import collections
import csv

import _ctx  # noqa: F401

from plasmidann.cascade import (
    check_thresholds,
    classify,
    completeness,
    dark_evidence,
    explained_fraction,
    is_informative,
    n_dark_databases,
    uninformative_spans,
)
from plasmidann.plasmidscope import TIER as PS_TIER
from plasmidann.plasmidscope import annot_label

cfg = snakemake.params.thresholds
# Fail loudly and early on an incoherent threshold block rather than producing a table
# nobody can interpret. The JSON schema checks ranges; this checks relationships.
check_thresholds(cfg)
tier_order = snakemake.params.tier_order

# ------------------------------------------------------------------------------------
# Gather every hit for every protein, from every tier.
# ------------------------------------------------------------------------------------
by_query = collections.defaultdict(list)
unnamed = collections.defaultdict(list)
named = collections.defaultdict(list)
for f in snakemake.input.hits:
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["coverage"] = float(r["coverage"]) if r["coverage"] else ""
            r["target_coverage"] = float(r["target_coverage"]) if r["target_coverage"] else None
            # Empty for a tier that reports no alignment span (pharokka). Left as the
            # empty string rather than coerced: cascade.classify reads '' as "no span".
            r["start"] = int(r["start"]) if r["start"] else ""
            r["end"] = int(r["end"]) if r["end"] else ""
            by_query[r["query"]].append(r)
            if not is_informative(r["label"]):
                unnamed[r["query"]].append(r)
            else:
                named[r["query"]].append(r)

# The informative explained fraction is taken from the final tier's cumulative spans,
# which were merged as the cascade descended. qlen comes with it. Each tier writes forward
# everything it inherited, so the last tier's spans hold every protein the cascade ever
# saw.
explained, qlen = {}, {}
with open(snakemake.input.spans, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        explained[r["seq_id"]] = float(r["explained_fraction"])
        qlen[r["seq_id"]] = int(r["qlen"])

seq_ids = [l[1:].split()[0] for l in open(snakemake.input.faa) if l[0] == ">"]

with open(snakemake.input.selection, newline="") as fh:
    selection = list(csv.DictReader(fh, delimiter="\t"))
# Only representatives with members need their row kept for copying.
copied_from = {r["search_representative"] for r in selection if r["role"] == "member"}
# AntiFam-flagged proteins skip Tier 0 as well: PlasmidScope's name for one is not used.
antifam = {r["seq_id"] for r in selection if r["role"] == "artefact_antifam"}

cols = [
    "seq_id",
    # where the result comes from; see the module docstring
    "annot_source", "annot_representative",
    # what it is
    "annot_tier", "annot_label", "functional_class", "homology_depth",
    "annot_qcov", "annot_tcov", "annot_evalue", "n_informative_hits",
    # annot_label comes from the most authoritative tier that named the protein; the
    # label with the best E-value across tiers is kept beside it (cascade.classify).
    "best_evalue_label", "best_evalue_tier",
    # 1 when every informative name is domain-level ("X domain-containing protein"), which
    # makes the protein DOMAIN_ONLY however much of it those names cover.
    "named_by_domain_only",
    # EVERY informative label, not just the winning one. The backbone stop-list matches
    # Pfam family names, and only the Pfam tiers emit those; when labels were ranked by
    # E-value alone an nr hit with free text routinely took annot_label away from a curated
    # Pfam assignment on the same protein. A complete replication initiator carrying
    # RepA_N at T1 then reached the screening pool with the guard that exists to stop it
    # never firing, because the only column S5 could read no longer held a Pfam name.
    "informative_labels", "informative_tiers",
    # how much of it is accounted for
    "explained_fraction", "annot_completeness", "meets_min_explained",
    # what the dark evidence says
    "dark_covered_fraction", "dark_completeness", "dark_evidence",
    "n_dark_databases", "uninformative_labels", "uninformative_tiers",
    # provenance: every row carries the thresholds that produced it (P4)
    "thr_min_coverage", "thr_min_explained", "thr_narrow_at", "thr_full_at", "thr_partial_at",
]

with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    kept = {}
    for sid in seq_ids:
        ef = explained.get(sid, 0.0)
        length = qlen.get(sid, 0)
        u = unnamed.get(sid, [])
        n = named.get(sid, [])
        labels = [x["label"] for x in u]

        # Dark coverage: the same merge, over the spans of hits that named nothing.
        dcf = explained_fraction(length, uninformative_spans(
            [x for x in u if x["start"] != "" and x["end"] != ""]))

        classified = classify(by_query.get(sid, []), explained=ef,
                              min_coverage=cfg["min_coverage"], tier_order=tier_order)
        span_measured = classified.pop("span_measured")

        row = {
            "seq_id": sid,
            "annot_source": "self",
            **classified,
            "explained_fraction": ef,
            # NOT_MEASURED when the class rests on a family-level assignment alone: the
            # explained fraction is then 0 because nothing measured it, and reporting
            # completeness NONE would read as "nothing matched" on a protein whose whole
            # family is known. Spec section 2.9: absence of a measurement is a status.
            "annot_completeness": (
                completeness(ef, full_at=cfg["full_at"], partial_at=cfg["partial_at"])
                if span_measured else "NOT_MEASURED"),
            # min_explained as a reported flag rather than a filter: the protein is in the
            # table either way, and this column can be recomputed at any threshold.
            "meets_min_explained": int(ef >= cfg["min_explained"]),
            "dark_covered_fraction": dcf,
            "dark_completeness": completeness(dcf, full_at=cfg["full_at"],
                                              partial_at=cfg["partial_at"]),
            "dark_evidence": dark_evidence(labels),
            "n_dark_databases": n_dark_databases(u),
            # The labels themselves are kept: someone else has described this protein, and
            # that is evidence FOR screening it, not against.
            "informative_labels": " | ".join(dict.fromkeys(x["label"] for x in n)),
            "informative_tiers": ",".join(dict.fromkeys(x["tier"] for x in n)),
            "uninformative_labels": " | ".join(dict.fromkeys(labels)),
            "uninformative_tiers": ",".join(dict.fromkeys(x["tier"] for x in u)),
            "thr_min_coverage": cfg["min_coverage"],
            "thr_min_explained": cfg["min_explained"],
            "thr_narrow_at": cfg["narrow_at"],
            "thr_full_at": cfg["full_at"],
            "thr_partial_at": cfg["partial_at"],
        }
        w.writerow(row)
        if sid in copied_from:
            kept[sid] = row

    # Proteins PlasmidScope annotates never entered the cascade (prepare_control); they
    # are FUNCTIONAL on its eggNOG result. eggNOG reports no alignment span, so nothing
    # measured how much of the protein is explained: the span fields stay empty and
    # completeness is NOT_MEASURED, as for any family-level assignment.
    n_ps = 0
    with open(snakemake.input.ps, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["ps_class"] != "ANNOTATED" or r["seq_id"] in antifam:
                continue
            n_ps += 1
            w.writerow({
                "seq_id": r["seq_id"],
                "annot_source": "plasmidscope",
                "annot_tier": PS_TIER,
                "annot_label": annot_label(r),
                "functional_class": "FUNCTIONAL",
                "n_informative_hits": 1,
                "informative_labels": annot_label(r),
                "informative_tiers": PS_TIER,
                "annot_completeness": "NOT_MEASURED",
                "thr_min_coverage": cfg["min_coverage"],
                "thr_min_explained": cfg["min_explained"],
                "thr_narrow_at": cfg["narrow_at"],
                "thr_full_at": cfg["full_at"],
                "thr_partial_at": cfg["partial_at"],
            })

    # Members of a search cluster take their representative's row; the rest of the
    # unique proteins were never searched, and say so.
    n_member = n_not_searched = n_antifam = 0
    for r in selection:
        if r["role"] == "member":
            n_member += 1
            w.writerow({**kept[r["search_representative"]], "seq_id": r["seq_id"],
                        "annot_source": "representative",
                        "annot_representative": r["search_representative"]})
        elif r["role"] == "not_selected":
            n_not_searched += 1
            w.writerow({"seq_id": r["seq_id"], "annot_source": "not_searched",
                        "functional_class": "NOT_SEARCHED"})
        elif r["role"] == "artefact_antifam":
            n_antifam += 1
            w.writerow({"seq_id": r["seq_id"], "annot_source": "artefact_antifam",
                        "functional_class": "NOT_SEARCHED"})

print(f"resolved {len(seq_ids)} proteins across tiers {tier_order}; "
      f"{n_member} members took their representative's result; "
      f"{n_ps} more resolved by PlasmidScope; {n_not_searched} not searched; "
      f"{n_antifam} AntiFam-flagged, not searched")
