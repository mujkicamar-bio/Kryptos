"""S0b: the clonal registry, one row per analysis-set plasmid.

Records per plasmid its MOB-suite cluster (a plasmid without one is its own singleton
cluster), the observed host and the source that named it, lifestyle (isolate or
metagenome, so a count can be restricted to either), MOB-suite's predicted host range,
topology, size and habitat. Independence is not counted here: the downstream counts use the
Mash lineages of plasmid_lineage.tsv, and the MOB cluster gives only the MOB_count and the
SINGLE_MOB and CROSS_MOB labels.

The host comes from three sources in order of preference (plasmidann.hosts): PLSDB species,
PlasmidScope's per-record host (which carries IMG/PR's) and the GenBank/RefSeq source
organism. Together they name a host for 61.0% of the analysis set (98.4% of isolate
plasmids, 23.9% of metagenomic ones), against 31.8% for PLSDB alone.

predicted_host_range is MOB-suite's predicted host range (master table mob_host_range): the
taxa in which plasmids with the same replicon and relaxase types have been seen, at
whatever rank fits (a genus, a family, an order, one or several phyla). It is a prediction
and a range, not an observed host, so it has its own column and never enters species,
genus or the host counts. For the 55,936 plasmids with no observed host it gives a range
for 18,244, almost all metagenomic.
"""
import csv

import _ctx  # noqa: F401

from plasmidann import hosts

keep = {l.strip() for l in open(snakemake.input.ids) if l.strip()}

ps_host = {}
with open(snakemake.input.ps_hosts, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["plasmid_id"] in keep:
            ps_host[r["plasmid_id"]] = r.get("host", "")
organism, lifestyle = {}, {}
with open(snakemake.input.working_set, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["plasmid_id"] in keep:
            organism[r["plasmid_id"]] = r.get("organism", "")
            lifestyle[r["plasmid_id"]] = r.get("lifestyle", "")

n, n_missing, n_host = 0, 0, 0
with open(snakemake.input.master, newline="") as fh, \
        open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=["plasmid_id", "mob_cluster", "species", "genus",
                                        "host_source", "predicted_host_range", "lifestyle",
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
        pid = r["plasmid_id"]
        species, genus, source = hosts.resolve([
            ("plsdb", r.get("plsdb_species", "")), ("plasmidscope", ps_host.get(pid, "")),
            ("organism", organism.get(pid, ""))])
        n_host += bool(genus)
        w.writerow({"plasmid_id": pid, "mob_cluster": cluster, "species": species,
                    "genus": genus, "host_source": source,
                    "predicted_host_range": (r.get("mob_host_range") or "").strip(),
                    "lifestyle": lifestyle.get(pid, ""), "topology": r.get("topology", ""),
                    "size_bp": r.get("size_bp", ""), "hab_top": r.get("hab_top", "")})
        n += 1

print(f"registry: {n} plasmids, {n_missing} without a MOB cluster (own lineage), "
      f"{n_host} with a host")
assert n == len(keep), f"registry has {n} rows for {len(keep)} analysis-set plasmids"
