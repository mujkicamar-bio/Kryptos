"""Build a DRAFT category configuration from the observed label vocabulary.

WHEN THIS RUNS

After a full annotation run, once results/s4c/protein_labels.tsv exists. It is not part of
the pipeline DAG: it is a one-off that produces something for a person to read.

WHAT IT DOES, AND WHY IT STOPS WHERE IT DOES

It applies the patterns in config/category_mining.yaml to the DESCRIPTIONS of the observed
labels and writes two files:

  category_draft.yaml    every candidate member of every category, with the number of
                         proteins carrying it, in the shape config/label_categories.yaml
                         expects
  label_frequency.tsv    the full observed vocabulary, per kind, ranked by protein count

It does NOT write config/label_categories.yaml. The draft is reviewed first, and the
reviewed result is what ships. Measured on the installed Pfam-A 38.2, the reason for that
separation:

  /replicat/    matches  67 family descriptions   (the deleted hand list named 16)
  /conjug/      matches  42                        (it named 15)
  /transposas/  matches  76                        (it named 15)
  /toxin/       matches 317  - including ABC_toxin_N, an insect toxin, and ADPRTs_Tse2, a
                               T6SS effector; neither is a toxin-antitoxin system

So the pattern reaches a scope no hand list can, and the review removes what the pattern
cannot judge. Applying the pattern at query time instead would keep the scope and lose the
review - and a substring match against a live vocabulary is how a mislabelled neighbourhood
sends someone to the bench to test the wrong thing.

The protein count beside each candidate is the evidence for the review: it says what
keeping or dropping that candidate costs, before the decision is made.
"""
import _ctx  # noqa: F401
import collections
import csv
import re

import yaml

from plasmidann import normalise, pfam_meta

# Label kinds whose text is free-form and must be normalised before it is counted or
# matched. The controlled kinds are exact strings and are used as written.
_FREE_TEXT = {"pgap_product", "swissprot_product", "eggnog_description"}

mining = (yaml.safe_load(open(snakemake.input.mining).read()) or {}).get("mine") or {}
pfam = pfam_meta.load(snakemake.input.pfam_dat)

# ------------------------------------------------------------------------------------
# The observed vocabulary: how many distinct proteins carry each (kind, label).
#
# Proteins, not rows: a label seen on one protein through three tiers is one protein's
# worth of evidence, and counting rows would make widely searched labels look commoner.
# ------------------------------------------------------------------------------------
proteins = collections.defaultdict(set)
with open(snakemake.input.labels, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        label = r["label"]
        if r["kind"] in _FREE_TEXT:
            label = normalise.product_name(label)
        if label:
            proteins[(r["kind"], label)].add(r["protein_id"])


def describe(kind, label):
    """The text a pattern is matched against for one label.

    For a Pfam family this is the release description, which is where the functional words
    live: 'RepA_N' contains no matchable word, 'Replication initiator protein A (RepA)
    N-terminus' does. For every other kind the label is its own description.
    """
    if kind == "pfam_family":
        return pfam.get(label, {}).get("description", "") or label
    return label


with open(snakemake.output.frequency, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["kind", "label", "n_proteins", "description"])
    for (kind, label), ids in sorted(proteins.items(),
                                     key=lambda kv: (-len(kv[1]), kv[0])):
        description = describe(kind, label) if kind == "pfam_family" else ""
        w.writerow([kind, label, len(ids), description])

# ------------------------------------------------------------------------------------
# Candidates per category.
# ------------------------------------------------------------------------------------
draft = {}
for category, spec in mining.items():
    spec = spec or {}
    description_pattern = re.compile(spec["description"], re.I) \
        if spec.get("description") else None
    symbol_pattern = re.compile(spec["symbol"]) if spec.get("symbol") else None
    cog = set(spec.get("cog_category") or [])

    by_kind = collections.defaultdict(list)
    for (kind, label), ids in proteins.items():
        if kind == "cog_category":
            matched = label in cog
        elif kind == "gene_symbol":
            matched = bool(symbol_pattern and symbol_pattern.search(label))
        else:
            matched = bool(description_pattern
                           and description_pattern.search(describe(kind, label)))
        if matched:
            by_kind[kind].append((len(ids), label))

    if by_kind:
        draft[category] = {
            "exact": {kind: [label for _, label in sorted(values, reverse=True)]
                      for kind, values in sorted(by_kind.items())},
            "_candidate_counts": {
                kind: {label: n for n, label in sorted(values, reverse=True)}
                for kind, values in sorted(by_kind.items())},
        }

header = (
    "# DRAFT - generated by workflow/scripts/build_categories.py. NOT the live config.\n"
    "#\n"
    "# Every entry below is a CANDIDATE produced by a pattern in\n"
    "# config/category_mining.yaml. Patterns reach a scope no hand-written list can, and\n"
    "# they also match things that do not belong: measured on Pfam-A 38.2, /toxin/ hits\n"
    "# 317 family descriptions of which ABC_toxin_N is an insect toxin and ADPRTs_Tse2 is\n"
    "# a T6SS effector.\n"
    "#\n"
    "# Review this file, delete what does not belong, delete the _candidate_counts blocks,\n"
    "# and copy the result into config/label_categories.yaml. Record for each category, in\n"
    "# docs/PARAMETER_PROVENANCE.md: the pattern that produced it, what was removed and\n"
    "# why, and the number of proteins it ends up covering.\n"
    "#\n"
    "# _candidate_counts gives the number of distinct proteins carrying each candidate, so\n"
    "# the cost of keeping or dropping one is visible before the decision is made.\n")

with open(snakemake.output.draft, "w") as out:
    out.write(header)
    yaml.safe_dump({"categories": draft}, out, default_flow_style=False, sort_keys=True)

n_proteins = len({p for ids in proteins.values() for p in ids})
print(f"build_categories: {len(proteins)} distinct (kind, label) pairs over "
      f"{n_proteins} proteins")
for category in sorted(draft):
    counts = draft[category]["_candidate_counts"]
    covered = len({p for kind, labels in counts.items()
                   for label in labels for p in proteins[(kind, label)]})
    per_kind = ", ".join(f"{k}={len(v)}" for k, v in sorted(counts.items()))
    print(f"  {category:<26} {sum(len(v) for v in counts.values()):>5} candidates "
          f"({per_kind}), covering {covered} proteins")
