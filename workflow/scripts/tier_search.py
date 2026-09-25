"""S3: search one tier, then hand the next tier only what stayed unexplained.

This script is where a search tool's raw output becomes evidence. Four things happen, in
this order, and the order matters:

  1. Run the tier's tool over whatever the previous tier could not explain.
  2. Reject hits that are not significant, BEFORE they can contribute anything.
  3. Accumulate INFORMATIVE spans (only) into the running explained fraction.
  4. Write out the residue for the next tier, honouring the sweep cohort.

Each step carries the reasoning for why it is done the way it is, including the measured
consequence of the alternative where one was measured.
"""
import csv
import os
import pathlib
import subprocess

import _ctx  # noqa: F401  - puts src/ on sys.path for the plasmidann package

from plasmidann import pharokka, scratch
from plasmidann.cascade import (
    as_float,
    explained_fraction,
    is_informative,
    narrow_by_explained,
    passes_significance,
)

spec = snakemake.params.spec              # one entry from cascade.yaml: tiers
faa = snakemake.input.faa                 # this tier's queries
max_evalue = spec["max_evalue"]           # None for bit-score tiers such as T1 (--cut_ga)
hmmer_z = snakemake.params.hmmer_z        # fixed -Z / --domZ, see docs/annotation_statistics.md
max_target_seqs = snakemake.params.max_target_seqs

ids = [l[1:].split()[0] for l in open(faa) if l[0] == ">"]

# Proteins an earlier tier already NAMED from a source this tier defers to are not searched
# here (skip_if_named_by in config/cascade.yaml). Decided for nr: a protein Pfam or
# Swiss-Prot named keeps that curated name, and nr - free text, one representative per
# cluster - is spent only on what they could not name. They still pass through the
# narrowing below with the explained fraction they arrived with; the sweep cohort is not
# exempt, since the rule is about which database may name a protein, not about narrow_at.
skip_sources = set(spec.get("skip_if_named_by") or [])
skipped = set()
if skip_sources:
    id_set0 = set(ids)
    for path in snakemake.input.named:
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["informative"] == "True" and r["source"] in skip_sources \
                        and r["query"] in id_set0:
                    skipped.add(r["query"])
search_ids = [i for i in ids if i not in skipped]

# Proteins that bypass narrowing entirely and are searched by EVERY tier. For these the
# counterfactual actually exists, which is the only way to measure what the narrowing
# threshold cost. Empty file (or fraction 0) disables the cohort.
sweep_ids = set()
if getattr(snakemake.input, "sweep", None):
    sweep_ids = {l.strip() for l in open(snakemake.input.sweep) if l.strip()}

# ------------------------------------------------------------------------------------
# Carry forward what earlier tiers already explained.
#
# spans hold ONLY informative alignments: a 'hypothetical protein' hit explains nothing,
# however well it aligns, and must never satisfy the narrowing threshold - another
# database may still name the protein.
# ------------------------------------------------------------------------------------
spans, qlen = {}, {}
prev = snakemake.input.spans
for f in ([prev] if isinstance(prev, str) else list(prev)):
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            qlen[r["seq_id"]] = int(r["qlen"])
            spans[r["seq_id"]] = [tuple(map(int, iv.split("-")))
                                  for iv in filter(None, r["intervals"].split(";"))]

named = {}     # EVERY informative hit per query, this tier
best = {}      # which of them supplies the winning label
unnamed = {}   # every uninformative hit, kept as evidence someone else has seen this protein
n_rejected = 0 # hits dropped for insignificance - reported, so the gate is visible


def record(q, label, qcov, tcov, ev, start, end, tlen, accession="", source=None,
           category="", identity="", align_length="", bitscore="", target_length=""):
    """Admit one parsed hit, or reject it.

    `start`/`end` may be None for a tier that reports no alignment span - the pharokka
    tier, whose raw alignments the tool deletes on exit. Such a hit contributes no span
    to `explained`; it is a family-level assignment, and cascade.classify treats it as one.
    The columns are written empty rather than invented.

    SIGNIFICANCE GATE. Nothing enters until it clears the tier's declared max_evalue.
    This matters because an insignificant domain that contributed a span could push a
    protein over the narrowing threshold and so withhold it from every deeper tier.
    Measured: 7.2% of T2 resolutions depend on domains that are not individually
    significant, because -E sets only the SEQUENCE threshold.

    identity, align_length, bitscore and target_length are what a tool reports and are
    written empty where it reports none (spec sections 19 and 21 require them for the
    DIAMOND tiers). Identity is what any later limit on label transfer needs, and it cannot
    be recovered from hits.tsv without repeating the search.
    """
    global n_rejected
    if not passes_significance(ev, max_evalue):
        n_rejected += 1
        return

    qlen[q] = tlen
    spanned = start is not None and end is not None
    if is_informative(label):
        if spanned:
            spans.setdefault(q, []).append((start, end))
        # EVERY informative hit is kept, not one per tier. Two documented claims depend
        # on this, and neither holds if the remaining hits are discarded:
        #
        #   * n_informative_hits is meant to tell a reader whether 0.9 coverage came from
        #     one domain or six. Capped at one per tier, its maximum was the tier count.
        #   * the coordinates let explained_fraction be recomputed from hits.tsv as an
        #     independent cross-check of spans.tsv. Keeping only one hit per tier makes
        #     that recomputation smaller every time, so the check fails on exactly the
        #     multi-domain proteins it exists for - a replication initiator carrying
        #     RepA_N and Bac_RepA_C being the case that matters most here.
        hit = {"query": q, "label": label, "target_accession": accession,
               "coverage": round(qcov, 4) if qcov is not None else "",
               "target_coverage": round(tcov, 4) if tcov is not None else "",
               "evalue": ev, "informative": True,
               "start": start if spanned else "", "end": end if spanned else "",
               "source": source or spec["source"], "category": category,
               "identity": identity, "align_length": align_length, "bitscore": bitscore,
               "target_length": target_length}
        named.setdefault(q, []).append(hit)
        # The best hit is the most SIGNIFICANT one, not the widest: a longer alignment is
        # not a better identification. Measured: 15.4% of labels change. The case that
        # settled it was
        # ABC_membrane at E=1e-23 being chosen over Peptidase_C39 at E=6.5e-40 purely
        # because it aligned further.
        if q not in best or as_float(ev) < as_float(best[q]["evalue"]):
            best[q] = hit
    else:
        # Uninformative hits keep their COORDINATES now, not just their label. Those
        # coordinates become dark_covered_fraction at cascade_resolve: 95% of a protein
        # covered by 'hypothetical protein' across three databases is a real, conserved,
        # full-length unnamed protein; one 20-aa fragment hit is not. Discarding the spans
        # makes those two indistinguishable.
        unnamed.setdefault(q, []).append(
            {"query": q, "label": label, "target_accession": accession,
             "coverage": round(qcov, 4) if qcov is not None else "",
             "target_coverage": round(tcov, 4) if tcov is not None else "",
             "evalue": ev, "informative": False,
             "start": start if spanned else "", "end": end if spanned else "",
             "source": source or spec["source"], "category": category,
             "identity": identity, "align_length": align_length, "bitscore": bitscore,
             "target_length": target_length})


if search_ids:
    # A stable name: each tier has its own output directory, so no two jobs share it, and
    # a killed attempt's directory (a partial res.m8 of the nr tier is ~GBs) is cleared by
    # the rerun instead of accumulating beside the new one.
    tmp = scratch.scratch_dir(pathlib.Path(snakemake.output.hits).parent, "search_tmp")
    query = faa
    if skipped:
        query = f"{tmp}/query.faa"
        with open(query, "w") as out:
            emit = False
            for line in open(faa):
                if line[0] == ">":
                    emit = line[1:].split()[0] not in skipped
                if emit:
                    out.write(line)

    if spec["method"] == "hmmer":
        raw = f"{tmp}/dom.tbl"
        # -Z / --domZ pin the search space. hmmsearch reports
        # E = (sequences searched) x P(score | null), and by default the first term is the
        # actual input size. Since each tier's input is the previous tier's residue, an
        # unpinned -E would mean a different significance on
        # every tier. Pinning it makes an E-value mean "expected false positives across
        # the whole study". See docs/annotation_statistics.md section 3.
        cmd = (f"hmmsearch {spec['args']} -Z {hmmer_z} --domZ {hmmer_z} "
               f"--noali --cpu {snakemake.threads} --domtblout {raw} {spec['db']} {query}")
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)

        # --domtblout column map (0-indexed). Verified against real output in review;
        # every index below was independently confirmed correct.
        #   0  target name  = the PROTEIN (hmmsearch searches profiles against sequences)
        #   2  tlen         = protein length
        #   3  query name   = the Pfam family
        #   5  qlen         = profile length
        #   6  full-sequence E-value
        #  12  i-Evalue     = independent E-value for THIS domain  <- the governing statistic
        #  13  domain score (bits)
        #  15,16 hmm from,to
        #  17,18 ali from,to = coordinates on the protein
        for line in open(raw):
            if line.startswith("#"):
                continue
            f = line.split()
            tlen, hlen = int(f[2]), int(f[5])
            a, b = int(f[17]), int(f[18])
            # The i-Evalue is used, not the sequence E-value. f[6] describes the whole
            # sequence;
            # f[12] is the statistic that governs an individual domain and is what HMMER's
            # own documentation directs users to trust when a sequence has several.
            # f[4] is the QUERY accession, which for hmmsearch is the Pfam accession
            # (PF06970.19) - hmmsearch searches profiles against sequences, so the "query"
            # is the profile. f[1] is the target accession, which is the protein's and is
            # '-' for our FASTA. The accession is kept because it is stable across Pfam
            # releases while a family NAME can be changed by a curator.
            record(f[0], f[3],
                   (b - a + 1) / tlen if tlen else 0.0,
                   (int(f[16]) - int(f[15]) + 1) / hlen if hlen else 0.0,
                   f[12], a, b, tlen, accession=f[4], bitscore=f[13],
                   target_length=hlen)

    elif spec["method"] == "pharokka":
        # The phage tier (spec section 18): pharokka in protein mode, used as its authors
        # ship it. It runs in its own environment, named by `exe`; pharokka checks for its
        # helpers on PATH, so that environment's bin directory is put there for the call.
        # See plasmidann.pharokka for what the output does and does not carry.
        exe = pathlib.Path(spec["exe"]).resolve()
        outdir = f"{tmp}/pharokka"
        evalue_flag = "" if max_evalue is None else f"-e {max_evalue} "
        subprocess.run(
            f"{exe} proteins -i {query} -o {outdir} -d {spec['db']} -t {snakemake.threads} "
            f"-f -p tier {evalue_flag}{spec.get('args', '')}",
            shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            env={**os.environ, "PATH": f"{exe.parent}:{os.environ.get('PATH', '')}"})

        merged = pathlib.Path(outdir) / "tier_full_merged_output.tsv"
        if not merged.exists():
            raise SystemExit(
                f"[{spec['id']}] pharokka exited 0 but wrote no {merged.name}. See "
                f"{outdir}/logs; the scratch directory is kept.")
        for hit in pharokka.parse_merged(merged.read_text()):
            # No span: start and end are None, and record() writes them empty. The label
            # goes through is_informative exactly as at every other tier, so a family
            # named 'hypothetical protein' does not make a dark protein annotated.
            record(hit["query"], hit["label"], None, None, hit["evalue"], None, None,
                   hit["query_length"], accession=hit["family_id"],
                   source=hit["source"], category=hit["category"])

    elif spec["method"] == "diamond":
        raw = f"{tmp}/res.m8"
        # The DIAMOND threshold is DECLARED, never left to the tool. DIAMOND's own default
        # is 0.001, which across 3.5M queries against nr - at the deepest tier, where a
        # spurious hit permanently removes a genuine dark protein from the pool - is a very
        # large number of expected false positives.
        #
        # The flag is OMITTED rather than stringified when a tier declares max_evalue: null,
        # meaning "the tool's own threshold decides", the convention T1 uses. Interpolating
        # None unconditionally emits `--evalue None`, which DIAMOND does not reject: it
        # parses to 0, returns ZERO hits and exits 0. On the nr tier that reads as "nothing
        # in nr matched any of 3.5 million proteins", and every downstream stage accepts it.
        # Verified: `--evalue None` gives 0 rows and exit 0 against a database containing
        # the query itself.
        # sseqid is the subject ACCESSION and stitle is free text. Both are kept: the
        # title is the only human-readable identification at T4, and the accession is the
        # only key that joins to UniProt keywords or to RefSeq. Measured on the NCBI
        # swissprot database in use here, the title reads 'P62554.1 RecName: Full=Toxin
        # CcdB; ... [Escherichia coli]' - NCBI's rendering, with no gene symbol - so the
        # accession carries the whole join.
        evalue_flag = "" if max_evalue is None else f"--evalue {max_evalue} "
        # --tmpdir sends DIAMOND's spill space to node-local disk. Left unset, DIAMOND
        # defaults it to the output file's directory, which is on the shared filesystem:
        # the nr tier then writes its intermediates across the network while streaming
        # 357 GB of database in the other direction. Measured on the 3,910-query benchmark,
        # the shared-filesystem run came in 15% over the local-disk prediction, and the
        # gap grows with the query count.
        cmd = (f"diamond blastp -q {query} -d {spec['db']} -o {raw} {spec['args']} "
               f"{evalue_flag}--threads {snakemake.threads} "
               f"--tmpdir {snakemake.resources.tmpdir} "
               f"--max-target-seqs {max_target_seqs} "
               f"--outfmt 6 qseqid sseqid stitle qcovhsp scovhsp evalue qstart qend qlen "
               f"pident length bitscore slen --quiet")
        subprocess.run(cmd, shell=True, check=True)
        for line in open(raw):
            (q, sid, title, qc, tc, ev, qs, qe, ql,
             pid, alen, bits, slen) = line.rstrip("\n").split("\t")
            record(q, title, float(qc) / 100, float(tc) / 100, ev,
                   int(qs), int(qe), int(ql), accession=sid, identity=float(pid) / 100,
                   align_length=int(alen), bitscore=bits, target_length=int(slen))

    else:
        # Previously `else` WAS the diamond branch, so a tier declaring any unrecognised
        # method ran DIAMOND against a database meant for another tool. That fails as a
        # database error at best and as a silently empty tier at worst.
        raise SystemExit(
            f"tier {spec['id']} declares method {spec['method']!r}, which is not one of "
            "hmmer, pharokka or diamond. A tier whose method is not recognised would "
            "otherwise fall through to another tool.")

# ------------------------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------------------------
# start/end are new columns: cascade_resolve needs the coordinates of UNINFORMATIVE hits
# to compute dark_covered_fraction, and they let explained_fraction be recomputed from
# hits.tsv as an independent cross-check.
# `is_best` marks the hit that supplies the winning label for this protein at this tier.
# It is a flag rather than a separate file because the other hits are not runners-up to be
# discarded: they are the rest of the protein's domain architecture, and the merged coverage
# that decides FUNCTIONAL versus DOMAIN_ONLY is built from all of them.
# `source` is the database the tier searched, declared per tier in config: it is what
# decides a label's kind downstream, and reading it from the row replaces guessing it from
# the tier's position. `category` is the functional group a database assigns to a family
# where it has one (pharokka, CARD); empty elsewhere.
cols = ["query", "label", "target_accession", "coverage", "target_coverage", "evalue",
        "informative", "is_best", "start", "end", "tier", "source", "category",
        "threshold", "max_evalue", "identity", "align_length", "bitscore", "target_length"]
with open(snakemake.output.hits, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    all_hits = ([x for v in named.values() for x in v]
                + [x for v in unnamed.values() for x in v])
    winners = {id(h) for h in best.values()}
    for h in all_hits:
        # Every row carries the thresholds that produced it (design principle P4), so any
        # downstream table can be traced back to the numbers that made it.
        w.writerow({**h, "is_best": int(id(h) in winners), "tier": spec["id"],
                    "threshold": spec["args"], "max_evalue": max_evalue})

explained = {q: explained_fraction(qlen.get(q, 0), spans.get(q, [])) for q in qlen}
with open(snakemake.output.spans, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["seq_id", "qlen", "intervals", "explained_fraction"])
    for q in qlen:
        w.writerow([q, qlen[q], ";".join(f"{a}-{b}" for a, b in spans.get(q, [])),
                    explained.get(q, 0.0)])

# The id set is built ONCE, outside the comprehension. Building it inside reconstructs it
# per explained id, which is quadratic: measured at 13.0 s for n=20,000 and projecting to
# 1.6-4.6 days at n=3.49M - spent after the search has already finished.
id_set = set(ids)
explained_here = {k: v for k, v in explained.items() if k in id_set}

# narrow_at, NOT min_explained. Narrowing decides what to keep SEARCHING; min_explained is
# applied post hoc at cascade_resolve and stays sweepable because of that separation.
keep = set(narrow_by_explained(ids, explained_here, snakemake.params.narrow_at))

# The sweep cohort is never narrowed away, whatever its explained fraction.
keep |= (sweep_ids & id_set)

with open(snakemake.output.unresolved, "w") as out:
    emit = False
    for line in open(faa):
        if line[0] == ">":
            emit = line[1:].split()[0] in keep
        if emit:
            out.write(line)

print(f"[{spec['id']}] received={len(ids)} skipped_named_by_{'_'.join(sorted(skip_sources)) or 'none'}"
      f"={len(skipped)} queried={len(search_ids)} informative={len(best)} "
      f"informative_hits={sum(len(v) for v in named.values())} unnamed={len(unnamed)} "
      f"rejected_insignificant={n_rejected} carried_forward={len(keep)} "
      f"(sweep_cohort={len(sweep_ids & id_set)})")

# The scratch directory is removed only here, on the ordinary path. A script that raised
# never reaches this line, and its intermediates - the raw domtbl or m8 the tool wrote -
# are what a tier failure is diagnosed from. Guarded: an empty tier created no directory.
if search_ids:
    scratch.release(tmp)
