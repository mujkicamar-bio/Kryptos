"""Status values shared across stages: the companion <field>_status of a nullable value
says why the value is missing ("not run", "too few lineages" ...), which a bare null
cannot. A status written by one stage only (family_evolution's dnds_status values) is
defined where it is written."""

NOT_RUN = "NOT_RUN"                  # the stage was not executed for this record
NO_HIT = "NO_HIT"                    # the search ran and returned nothing above threshold
TOO_FEW_LINEAGES = "TOO_FEW_LINEAGES"  # fewer independent lineages than the estimator needs
NO_CONTEXT = "NO_CONTEXT"            # no occurrence has a named neighbour to measure synteny on
SUCCESS = "SUCCESS"                  # a value is present and usable
NOT_MEASURED = "NOT_MEASURED"        # no input record carried what the value is counted from
