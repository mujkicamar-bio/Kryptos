"""S0b: the clonal registry - which plasmids are independent observations.

Every later count of "independent occurrences" is computed over MOB clusters, not raw
plasmids. Without this, a family found on forty plasmids may be one clone sequenced forty
times, and the multi-lineage evidence at S7 and S9 means nothing.

mob_cluster is populated for essentially every plasmid in the collection (208,245 of the
master table), whereas host species is available for only 45,743 (~22%) - which is why the
plasmid lineage, not the host, is the unit of independence here. It is also the better
criterion: a protein family found on two unrelated plasmid backbones is stronger evidence
of a mobile functional unit than one found in two host species that happen to share the
same plasmid.
"""
import _ctx  # noqa: F401
import csv

keep = {l.strip() for l in open(snakemake.input.ids) if l.strip()}

n, n_missing = 0, 0
with open(snakemake.input.master, newline="") as fh, \
        open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=["plasmid_id", "mob_cluster", "species",
                                        "topology", "size_bp", "hab_top"],
                       delimiter="\t")
    w.writeheader()
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["plasmid_id"] not in keep:
            continue
        cluster = (r.get("mob_cluster") or "").strip()
        if not cluster or cluster == "NA":
            # A plasmid with no MOB assignment is its own lineage. Treating it as unknown
            # and dropping it would bias breadth counts toward well-typed plasmids.
            cluster = f"singleton:{r['plasmid_id']}"
            n_missing += 1
        w.writerow({"plasmid_id": r["plasmid_id"], "mob_cluster": cluster,
                    "species": r.get("plsdb_species", ""), "topology": r.get("topology", ""),
                    "size_bp": r.get("size_bp", ""), "hab_top": r.get("hab_top", "")})
        n += 1

print(f"registry: {n} plasmids, {n_missing} without a MOB cluster (own lineage)")
assert n == len(keep), f"registry has {n} rows for {len(keep)} analysis-set plasmids"
