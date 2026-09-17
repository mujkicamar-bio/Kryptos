"""Grouping labels into functional categories - the seam, not the grouping.

The grouping itself (replication, mobilisation, conjugation, ...) is derived from the label
vocabulary after a full annotation run, because that is the only way it can cover the
scope: Pfam-A 38.2 holds 30,134 families, of which 67 mention replication in their
description and 42 mention conjugation, where the hand list this replaces named 16 and 15.

What is built here is the mechanism and its default. With no rules configured, each
(kind, label) is its own category, so the enrichment machinery downstream is exercised on
real data from the first run and the grouping becomes a configuration change rather than a
rewrite.
"""
import pytest

from plasmidann import categories


def test_with_no_rules_a_label_is_its_own_category_qualified_by_kind():
    """The identity default. The kind qualifies the category name because 'repA' the gene
    symbol and 'RepA_N' the Pfam family are different statements, and a category that
    merged them could not be audited."""
    rules = categories.load_rules(None)

    assert categories.assign("pfam_family", "RepA_N", rules) == \
        [("pfam_family:RepA_N", "")]
    assert categories.assign("gene_symbol", "repA", rules) == [("gene_symbol:repA", "")]


def test_the_identity_map_reports_itself_as_the_identity():
    """A run whose categories are still the identity has not had the grouping applied, and
    a report must be able to say so rather than presenting per-label enrichment as if it
    were per-category."""
    assert categories.CategoryMap(categories.load_rules(None)).is_identity is True


def test_an_exact_rule_assigns_a_category():
    rules = categories.load_rules_from_text("""
categories:
  replication:
    exact:
      pfam_family: [RepA_N, Rep_3]
      gene_symbol: [repA]
""")

    assert categories.assign("pfam_family", "RepA_N", rules) == [("replication", "")]
    assert categories.assign("gene_symbol", "repA", rules) == [("replication", "")]


def test_a_label_matching_no_rule_falls_back_to_itself():
    """A label with no category must stay visible. Dropping it would make the vocabulary
    shrink silently as rules are added, which is how the previous list's nine dead entries
    went unnoticed for a whole version."""
    rules = categories.load_rules_from_text("""
categories:
  replication:
    exact:
      pfam_family: [RepA_N]
""")

    assert categories.assign("pfam_family", "MobA_MobL", rules) == \
        [("pfam_family:MobA_MobL", "")]


def test_a_label_may_belong_to_several_categories():
    """A relaxase is both mobilisation and, on a conjugative plasmid, part of transfer. A
    map that forced one category per label would decide that by rule order rather than by
    biology."""
    rules = categories.load_rules_from_text("""
categories:
  mobilisation:
    exact:
      pfam_family: [MobA_MobL]
  transfer:
    exact:
      pfam_family: [MobA_MobL]
""")

    assert sorted(categories.assign("pfam_family", "MobA_MobL", rules)) == \
        [("mobilisation", ""), ("transfer", "")]


def test_a_rule_naming_an_undeclared_kind_is_an_error():
    """A rule on a kind nothing emits can never fire, and it would look like a category
    with no members rather than a typo. Nine such entries in the list this replaces had
    never once matched anything."""
    with pytest.raises(ValueError, match="unknown label kind"):
        categories.load_rules_from_text("""
categories:
  replication:
    exact:
      pfam_familly: [RepA_N]
""")


def test_matching_is_exact_and_case_sensitive():
    """Substring matching would be a trap: 'PIN' and 'rve' are Pfam families whose names
    appear inside many unrelated descriptions, and a mislabelled protein sends someone to
    the bench to test the wrong thing."""
    rules = categories.load_rules_from_text("""
categories:
  toxin_antitoxin:
    exact:
      pfam_family: [PIN]
""")

    assert categories.assign("pfam_family", "PIN", rules) == [("toxin_antitoxin", "")]
    assert categories.assign("pfam_family", "PIN_2", rules) == [("pfam_family:PIN_2", "")]
    assert categories.assign("pfam_family", "pin", rules) == [("pfam_family:pin", "")]


def test_the_map_memoises_repeated_lookups():
    """The table has one row per protein per label; the same (kind, label) is looked up
    millions of times. A per-call scan of every rule would be the same quadratic class
    already fixed twice in this pipeline, in tier_search and in context_features."""
    cmap = categories.CategoryMap(categories.load_rules_from_text("""
categories:
  replication:
    exact:
      pfam_family: [RepA_N]
"""))

    first = cmap.categories_for("pfam_family", "RepA_N")
    second = cmap.categories_for("pfam_family", "RepA_N")

    assert first == second == [("replication", "")]
    assert cmap.categories_for("pfam_family", "RepA_N") is second


def test_a_subcategory_is_carried_but_does_not_change_the_category():
    """Flat for the test, finer for the record. A subcategory lets 'replication/initiator'
    and 'replication/control' be told apart in the table without splitting the category the
    significance test runs on - more categories means more tests and weaker FDR power, and
    re-cutting a finer grouping must not require re-running the pipeline."""
    rules = categories.load_rules_from_text("""
categories:
  replication:
    exact:
      pfam_family: [RepA_N, Rop]
    sub:
      pfam_family:
        RepA_N: initiator
        Rop: control
""")

    assert categories.assign("pfam_family", "Rop", rules) == [("replication", "control")]
    assert categories.assign("pfam_family", "RepA_N", rules) == \
        [("replication", "initiator")]


def test_a_label_with_a_subcategory_but_no_exact_rule_is_not_categorised():
    """A sub entry refines a membership that exact declares; it never creates one. If it
    did, the sub block would become a second, undocumented place where a category gains a
    member, and a reviewer reading the exact lists would not see the whole category."""
    rules = categories.load_rules_from_text("""
categories:
  replication:
    sub:
      pfam_family:
        RepA_N: initiator
""")

    assert categories.assign("pfam_family", "RepA_N", rules) == \
        [("pfam_family:RepA_N", "")]


def test_the_shipped_config_file_loads_and_is_still_the_identity():
    """config/label_categories.yaml ships with no rules, because the grouping is derived
    from a full annotation run. This test fails the day rules are added without updating
    it, which is the point at which a reviewer should look at them."""
    rules = categories.load_rules("config/label_categories.yaml")

    assert categories.CategoryMap(rules).is_identity is True, (
        "rules have been added to config/label_categories.yaml - update this test and "
        "record the provenance of each category in docs/PARAMETER_PROVENANCE.md")
