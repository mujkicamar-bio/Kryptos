"""The status vocabulary: the companion <field>_status of a nullable value says why the
value is missing ("not run", "too few members", "saturated" ...), which a bare null cannot."""

NOT_RUN = "NOT_RUN"                  # the stage was not executed for this record
NO_HIT = "NO_HIT"                    # the search ran and returned nothing above threshold
TOO_FEW_MEMBERS = "TOO_FEW_MEMBERS"  # below the configured minimum for the estimator
TOO_FEW_LINEAGES = "TOO_FEW_LINEAGES"  # fewer independent lineages than the estimator needs
NO_CONTEXT = "NO_CONTEXT"            # no occurrence has a named neighbour to measure synteny on
NO_DIVERGENCE = "NO_DIVERGENCE"      # sequences are identical; the statistic is undefined
SATURATED = "SATURATED"              # divergence too high for the estimate to be meaningful
NO_OUTPUT = "NO_OUTPUT"              # the tool exited successfully but produced nothing
SUCCESS = "SUCCESS"                  # a value is present and usable
NOT_MEASURED = "NOT_MEASURED"        # no input record carried what the value is counted from
