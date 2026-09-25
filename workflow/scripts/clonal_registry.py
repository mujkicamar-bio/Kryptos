"""S0b: the clonal registry - which plasmids are independent observations.

Every later count of "independent occurrences" is computed over MOB clusters, not raw
plasmids. Without this, a family found on forty plasmids may be one clone sequenced forty
times, and the multi-lineage evidence at S7 and S9 means nothing.

The plasmid lineage, not the host, is the unit of independence here: a protein family
found on two unrelated plasmid backbones is stronger evidence of a mobile functional unit
than one found in two host species that happen to share the same plasmid.

The host is recorded beside it, from three sources in order of preference - PLSDB species,
PlasmidScope's per-record host (which carries IMG/PR's), and the GenBank/RefSeq source
organism - which together name a host for 61.1% of the analysis set (98.6% of isolate
plasmids, 24.0% of metagenomic ones) against 31.9% for PLSDB alone (plasmidann.hosts). host_source says which source named it, and lifestyle whether
the plasmid came from an isolate or a metagenome, so a count can be restricted to either.

predicted_host_range is MOB-suite's predicted host range (master table mob_host_range):
the taxa in which plasmids with the same replicon and relaxase types have been seen, at
whatever rank fits - a genus, a family, an order, one or several phyla. It is a
prediction and a range, not an observed host, so it is kept in its own column and never
enters species, genus or the host counts (decided 2026-09-25). It is reported as a separate
measurement for every plasmid, hosted or not (recurrence, annotation_report). For the 55,785
plasmids with no observed host it adds a range for 18,170, almost all metagenomic.
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
