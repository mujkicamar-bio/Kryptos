"""The single status vocabulary (spec §7.2).

Every nullable numeric field in the pipeline carries a companion <field>_status drawn from
this set. The point is that "we did not run this", "this did not apply", "there were too
few sequences", "the estimate saturated" and "the tool crashed" are five different
scientific statements that a bare null would flatten into one.
"""

NOT_RUN = "NOT_RUN"                  # the stage was not executed for this record
NO_HIT = "NO_HIT"                    # the search ran and returned nothing above threshold
TOO_FEW_MEMBERS = "TOO_FEW_MEMBERS"  # below the configured minimum for the estimator
TOO_FEW_LINEAGES = "TOO_FEW_LINEAGES"  # fewer independent lineages than the estimator needs
NO_CONTEXT = "NO_CONTEXT"            # no occurrence has a named neighbour to measure synteny on
NO_DIVERGENCE = "NO_DIVERGENCE"      # sequences are identical; the statistic is undefined
SATURATED = "SATURATED"              # divergence too high for the estimate to be meaningful
NO_OUTPUT = "NO_OUTPUT"              # the tool exited successfully but produced nothing
FAILED = "FAILED"                    # the tool errored; see the log referenced in the row
NOT_APPLICABLE = "NOT_APPLICABLE"    # the field has no meaning for this record
SUCCESS = "SUCCESS"                  # a value is present and usable
NOT_MEASURED = "NOT_MEASURED"        # no input record carried what the value is counted from

ALL = frozenset({
    NOT_RUN, NO_HIT, TOO_FEW_MEMBERS, TOO_FEW_LINEAGES, NO_CONTEXT, NO_DIVERGENCE,
    SATURATED, NO_OUTPUT, FAILED, NOT_APPLICABLE, SUCCESS, NOT_MEASURED,
})
