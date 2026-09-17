"""Grouping tool labels into functional categories.

WHAT THIS IS AND IS NOT

This module is the mechanism for grouping labels, and its default is to do no grouping at
all. The biological categories - replication, mobilisation, conjugation, partition,
transposition, toxin-antitoxin - are NOT defined here. They are derived from the label
vocabulary after a full annotation run, from config/label_categories.yaml, and each one
needs a recorded justification like every other parameter in this pipeline.

The reason for that order is measured. The hand-written list this replaces named 16 of the
67 Pfam-A 38.2 families whose description mentions replication and 15 of the 42 that
mention conjugation, named no MobB and no MobD, and contained nine family names that do not
exist in Pfam-A at all - entries that had never matched anything in a whole version of the
pipeline, invisibly. Rules written against the labels the data actually contains can be
checked against the table; rules written from recollection cannot.

THE IDENTITY DEFAULT

With no rules configured, every (kind, label) is its own category, named 'kind:label'. Two
consequences, both wanted:

  * The enrichment machinery downstream runs on real data from the first annotation run,
    rather than waiting for a grouping that does not exist yet.
  * Applying the grouping later is a configuration change, not a rewrite.

The kind qualifies the category name because 'repA' as a gene symbol and 'RepA_N' as a Pfam
family are different statements about a protein. A category that merged them could not be
audited, and a reader could not tell which evidence produced it.

CATEGORY AND SUBCATEGORY

Every assignment is a (category, subcategory) pair. The CATEGORY is the flat axis the
significance test runs on; the SUBCATEGORY is carried beside it and never tested.

That split exists because the two things want opposite granularity. The test wants few
categories: every category tested is another hypothesis, and the false-discovery-rate
correction gets weaker as their number grows. The record wants many: 'replication' covering
both an initiator and a copy-number control protein loses a distinction someone reading the
table will want. Carrying the subcategory in the output means a finer grouping can be re-cut
from the table without re-running anything.

FALLBACK, NOT EXCLUSION

A label matching no rule falls back to being its own category. It is never dropped. A
vocabulary that shrank silently as rules were added is how the previous list's nine dead
entries survived unnoticed, and it is the one failure mode this design exists to prevent.
"""
import pathlib

import yaml

from plasmidann.labels import KINDS


def _check_kind(kind, category):
    """Reject a rule on a label kind nothing emits.

    Such a rule can never fire, and it would present as a category with no members rather
    than as the typo it is - which is exactly how nine entries of the list this replaces
    survived a whole version of the pipeline.
    """
    if kind not in KINDS:
        raise ValueError(
            f"unknown label kind {kind!r} in category {category!r}; "
            f"known kinds: {', '.join(sorted(KINDS))}")


def load_rules_from_text(text):
    """Parse category rules from YAML text.

    Shape:

        categories:
          replication:
            exact:
              pfam_family: [RepA_N, Rep_3]
              gene_symbol: [repA]
            sub:                      # optional
              pfam_family:
                RepA_N: initiator

    'sub' refines a membership that 'exact' declares and never creates one, so that there
    is exactly one place in the file where a category gains a member. A reviewer reading
    the exact lists sees the whole category.
    """
    document = yaml.safe_load(text) or {}
    declared = document.get("categories") or {}
    index = {}
    for category, spec in declared.items():
        spec = spec or {}
        subs = {}
        for kind, mapping in (spec.get("sub") or {}).items():
            _check_kind(kind, category)
            for value, name in (mapping or {}).items():
                subs[(kind, value)] = str(name)
        for kind, values in (spec.get("exact") or {}).items():
            _check_kind(kind, category)
            for value in values or []:
                index.setdefault((kind, value), []).append(
                    (category, subs.get((kind, value), "")))
    return {"exact": index}


def load_rules(path):
    """Load category rules from a path, or the empty (identity) rules when path is None."""
    if path is None:
        return {"exact": {}}
    return load_rules_from_text(pathlib.Path(path).read_text())


def assign(kind, label, rules):
    """The (category, subcategory) pairs of one label.

    Exact, case-sensitive matching. Substring matching would be a trap: 'PIN' and 'rve' are
    Pfam family names that appear inside many unrelated descriptions, and a mislabelled
    protein sends someone to the bench to test the wrong thing.

    A label may belong to several categories, and they are all returned. A relaxase is
    mobilisation and, on a conjugative plasmid, part of transfer; forcing one category per
    label would decide that by rule order rather than by biology.

    The subcategory is a finer name carried alongside the category and never tested. The
    significance test runs on the flat category, because every category tested is another
    hypothesis and the false-discovery-rate correction weakens as their number grows.
    """
    matched = rules["exact"].get((kind, label))
    return list(matched) if matched else [(f"{kind}:{label}", "")]


class CategoryMap:
    """A memoised view of the rules, for the per-row lookups.

    protein_labels has one row per protein per label, so the same (kind, label) is resolved
    millions of times. Scanning the rules per call would be the same quadratic class already
    fixed twice in this pipeline, in tier_search and in context_features.
    """

    def __init__(self, rules):
        self._rules = rules
        self._cache = {}

    @property
    def is_identity(self):
        """Whether no grouping is configured.

        A report has to be able to say that a run's categories are still per-label, rather
        than presenting per-label enrichment as though it were per-category.
        """
        return not self._rules["exact"]

    def categories_for(self, kind, label):
        key = (kind, label)
        if key not in self._cache:
            self._cache[key] = assign(kind, label, self._rules)
        return self._cache[key]
