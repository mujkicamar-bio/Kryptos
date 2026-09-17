"""S6b: cluster the dark set into families with MMseqs2.

A FAMILY IS A SEQUENCE CLUSTER. NOTHING ELSE.

An earlier draft of the design defined a family as ">=3 members from >=2 distinct MOB
clusters". That was wrong in the same way as baking a decision threshold into the search:
it makes the family COUNT a function of metadata completeness rather than of sequence, it
silently deletes lineage-restricted families when a novel system confined to one plasmid
lineage may be exactly what is worth finding, and it destroys the counterfactual - once a
cluster is dropped here, nothing downstream can ask what it would have looked like.

Size and distribution are ATTRIBUTES, computed here and weighted at S9 where they are
sweepable like every other weight. Clusters of one are retained and labelled ORPHAN.

30% identity / 50% coverage is FESNov's deep-homology setting, deliberately looser than the
70/80 used on 1.17 billion metagenomic sequences: we have 3.5 million, and at 70% real
families would fragment into singletons.
"""
import _ctx  # noqa: F401
import collections
import csv
import pathlib
import subprocess

cfg = snakemake.params.clustering
tmp = pathlib.Path(snakemake.output.tsv).parent / "mmseqs_tmp"
prefix = str(pathlib.Path(snakemake.output.tsv).parent / "dark")

subprocess.run(
    f"mmseqs easy-cluster {snakemake.input.faa} {prefix} {tmp} "
    f"--min-seq-id {cfg['min_seq_id']} -c {cfg['coverage']} "
    f"--cov-mode {cfg['cov_mode']} --cluster-mode {cfg['cluster_mode']} "
    f"--cluster-reassign --threads {snakemake.threads} -v 1",
    shell=True, check=True)

# easy-cluster writes <prefix>_cluster.tsv as representative<TAB>member.
members = collections.defaultdict(list)
with open(f"{prefix}_cluster.tsv") as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t")
        members[rep].append(mem)

# Plasmid provenance per protein, so family breadth can be measured over independent
# plasmid lineages rather than raw plasmid count. One clone sequenced forty times must
# count once, which is what the MOB-cluster registry is for.
seq_to_plasmids = collections.defaultdict(set)
seq_to_orf_count = collections.Counter()
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, ids = line.rstrip("\n").split("\t")
        oids = ids.split(",")
        seq_to_orf_count[sid] = len(oids)
        for oid in oids:
            seq_to_plasmids[sid].add(oid.rsplit("|", 1)[0])

mob_of = {}
with open(snakemake.input.registry, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        mob_of[r["plasmid_id"]] = r["mob_cluster"]

# n_members and n_orfs are DIFFERENT NUMBERS and both are needed.
#
# Clustering runs on the dereplicated set, so a protein whose sequence is identical on two
# hundred plasmids is ONE member: it clusters alone, is labelled ORPHAN, and fails the
# is_family reality test. That is a defensible definition - it is not a family of divergent
# homologs - but `n_members: 1` alone reads as "seen once", and a reader cannot then tell a
# genuine singleton from one of the most widely carried proteins in the collection.
# n_orfs is the number of actual gene copies behind the cluster.
cols = ["family_id", "representative", "n_members", "n_orfs", "n_plasmids",
        "n_mob_clusters", "family_class", "members"]
n_orphan = 0
with open(snakemake.output.tsv, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for i, (rep, mem) in enumerate(sorted(members.items()), start=1):
        plasmids = set()
        n_orfs = 0
        for m in mem:
            plasmids |= seq_to_plasmids.get(m, set())
            n_orfs += seq_to_orf_count.get(m, 0)
        mobs = {mob_of.get(p) for p in plasmids if mob_of.get(p)}
        # ORPHAN is a label, not an exclusion: a singleton is a legitimate category with a
        # known interpretation, and the orphan count is itself a diagnostic.
        family_class = "ORPHAN" if len(mem) == 1 else "FAMILY"
        n_orphan += family_class == "ORPHAN"
        w.writerow({"family_id": f"F{i:07d}", "representative": rep,
                    "n_members": len(mem), "n_orfs": n_orfs,
                    "n_plasmids": len(plasmids),
                    "n_mob_clusters": len(mobs), "family_class": family_class,
                    "members": ",".join(mem)})

print(f"families={len(members)} orphans={n_orphan} "
      f"multi_member={len(members) - n_orphan}")
