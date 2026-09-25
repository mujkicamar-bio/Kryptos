# Annotation statistics: thresholds, E-values and the two Pfam tiers

Companion to `plans/2026-09-10-pipeline-v2-design.md`. This document records *why* the search
parameters in `config/cascade.yaml` are what they are. Every number here is stamped into the
output rows it governs, so any table can be traced back to this file.

---

## 1. What `--cut_ga` is

Pfam curators set a **gathering threshold** (GA) for every family: a bit-score cutoff, chosen by
hand, such that sequences scoring above it are members of that family. The family's full
alignment is by definition the set of sequences above GA. `--cut_ga` tells `hmmsearch` to use
those per-family curated cutoffs instead of one E-value for everything.

**Why that is better than an E-value.** An E-value is a single number applied to every family
regardless of its length, its information content, or how prone it is to spurious matching. A
30-residue zinc finger and a 600-residue polymerase domain get the same bar. GA is a per-family
judgement by the person who curated that family — the difference between a speed limit set per
road and one blanket number for the whole country.

**What it buys.** Precision, and authority. A T1 hit is not our inference: it is Pfam's own
assertion that this protein belongs to this family.

**What it costs.** Recall, deliberately. GA is conservative, and the homologues it rejects are
concentrated among divergent sequences — which is most of the plasmid mobilome and nearly all
environmental sequence. For a dark-plasmidome project, that is precisely the population of
interest.

---

## 2. Why Pfam is searched twice

T1 and T2 use the same database and make **different claims with different authority**.

| tier | setting | claim |
|---|---|---|
| T1 | `--cut_ga` | *Pfam asserts this protein is a member of family X.* |
| T2 | `-E 1e-5 --domE 1e-5` | *This protein resembles family X below the curatorial bar.* |

Since tier order is authority order, keeping them separate means a T1 assignment can be reported
as a curated family membership, while a T2 assignment is clearly labelled as sub-GA inference
that can be down-weighted or stripped from any downstream analysis.

The alternatives are both worse. A single run at `-E 1e-5` collapses the two into an
undifferentiated pile in which a curated family assignment cannot be distinguished from a
marginal resemblance. A single run at `--cut_ga` discards the divergent homologues, which for
this project is discarding the subject matter.

**The two sets are not nested — and the direction is the opposite of the intuitive one.**
It is tempting to assume every T1 hit would also be found by T2, making T1 redundant. It is
false, and measuring how false it is changes what T2 can do.

Converting every Pfam-A gathering threshold to the E-value it corresponds to at
`-Z 3,497,616`, using each profile's own `STATS LOCAL FORWARD` parameters:

| | |
|---|---|
| Pfam-A profiles carrying GA and stats | 30,134 |
| profiles whose **GA is looser** than E=1e-5 at this Z | **29,074 (96.5%)** |
| median E-value at the curated GA | **7.6e-4** |
| median GA bit score | 27.0 |

So for 96.5% of Pfam families, `-E 1e-5` at this Z is **stricter** than the curator's own
cutoff. T1 finds hits T2 would reject, not the other way round.

### What that means for T2, and the decision taken

T2's stated job is to recover divergent homologues *below* the curatorial bar. It can only do
that for the 3.5% of families whose GA happens to be stricter than 1e-5. For the other 96.5%
there is no band between the two thresholds to search: anything T2 could find at E<=1e-5 has
already cleared GA and been caught by T1. So on this dataset T2 is the stricter of the two
Pfam tiers for most families, and much of what it returns is a re-find of T1's domains on
proteins that were not fully explained.

The cause is the size of the analysis set rather than anything about the tier. A gathering
threshold is a bit score a curator sets once, around 27 bits, and it does not move. An E-value
does: at `-Z 3,497,616` you need about 33 bits to reach 1e-5. Pinning `-Z` to the whole study -
which is correct, and is what makes E-values comparable across shards (section 3) - is what
lifted the bar above GA.

Pfam's second curated threshold does not help either. `--cut_nc`, the noise cutoff, is the
highest-scoring known non-member, so the NC-to-GA band is exactly the curated grey zone. All
30,134 profiles carry one, and the median GA minus NC gap is **0.2 bits**: Pfam sets GA
immediately above NC, so there is no usable band there.

**DECIDED: T2 stays at 1e-5.** Loosening it to ~1e-3 is the only setting that would open a real
sub-GA band, and its cost is precisely the error this project most needs to avoid - a false
positive at 3.5M queries permanently removes a genuine dark protein from the pool, with no
downstream stage able to recover it. Precision at the Pfam tiers is worth more here than sub-GA
recall. The measurement is kept because it is true and a reader should see it, not because the
setting is in question.

**Cost.** T2 searches only what T1 could not resolve, so it runs on a smaller input than T1 did.
Note that T2 will re-find the domains T1 already found on proteins that were not fully explained.
That is harmless — spans merge, and ties break to the shallowest tier — but it means T2's raw hit
count overstates what T2 contributed.

---

## 3. Why `-Z` is mandatory

### The mechanism

`hmmsearch` reports

```
E = (number of sequences searched) × P(score | null model)
```

and by default the first term is **the actual size of the input FASTA**. Measured on one
protein–profile pair, varying only the number of queries in the file:

| queries in input | reported E-value |
|---|---|
| 200 | 4.1e-36 |
| 2,000 | 4.1e-35 |
| 20,000 | 4.1e-34 |

A clean factor of ten per decade of input size, for an identical alignment.

### Why this breaks a cascade specifically

T2's input is *whatever T1 failed to resolve*. That quantity depends on the data in each shard
and on how T1 performed. Consequently:

1. **`-E 1e-5` means a different underlying significance on every shard.** Two shards with
   different T1 resolution rates get different effective stringency, so results are not
   comparable across shards of the same run.
2. **It changes whenever anything upstream changes.** Adjusting T1, the ORF length floor, or the
   dereplication step silently moves T2's threshold. The recorded number stays `1e-5` while the
   statistic it denotes moves.
3. **The threshold is therefore not a pre-registrable quantity**, which contradicts design
   principle P4.

### The fix

Pass a fixed, declared database size to every hmmer tier:

```
-Z 3497616 --domZ 3497616
```

3,497,616 is the number of unique protein sequences in the analysis set, from 9,317,050 ORFs
across 143,503 plasmids. With `-Z` fixed, an E-value means *the expected number of false
positives across the whole study* — a stable, interpretable, pre-registrable statement that does
not move when a shard boundary or an upstream stage changes. `--domZ` does the same for
domain-level E-values.

Without this, no E-value anywhere in the run means what it appears to mean.

---

## 4. Domain significance: `--domE` and the i-Evalue

A second, independent defect at T2.

`--cut_ga` applies gathering thresholds at **both** the sequence and the domain level, so weak
domains within a passing sequence are excluded by curation. `-E 1e-5` sets only the **sequence**
threshold. Every domain of a passing sequence then appears in `domtblout`, and the parser adds
all of them to the protein's explained span without checking whether any individual domain is
significant.

Consequence: a protein can be pushed over `min_explained` — and so withheld from every deeper
tier — by domains that are individually meaningless. Measured: **7.2% of T2 resolutions depend on
insignificant domains.**

There is also a statistic mismatch. The parser recorded field 7 of `domtblout`, the *full
sequence* E-value, when the quantity governing an individual domain is the **i-Evalue** at field
13 — the independent E-value, which accounts for multiple domains within one sequence and is the
value HMMER's own documentation directs users to trust.

### The fix

```
-E 1e-5 --domE 1e-5 --incdomE 1e-5 -Z 3497616 --domZ 3497616
```

plus an explicit i-Evalue filter applied at parse time, and recording the i-Evalue rather than
the sequence E-value for domain hits.

---

## 5. DIAMOND thresholds

DIAMOND's E-value threshold was never declared in v1. Its default is `--evalue 0.001`.

Across 3,497,616 queries against nr, that is a very large number of expected false positives, at
the deepest tier of the cascade — the last opportunity a protein has to remain dark. A spurious
hit there does maximum damage: it removes a genuine dark protein from the screening pool
permanently, and no downstream stage can recover it.

Every DIAMOND tier therefore declares `--evalue 1e-5` explicitly.

---

## 6. Where the E-value belongs in classification

Coverage answers *how much of this protein is named*. The E-value answers *how likely that name
is to be wrong*. Crossing them:

| | strong E-value | weak E-value |
|---|---|---|
| **high coverage** | `FUNCTIONAL` — correct | a spurious long alignment silently removes a real dark protein |
| **low coverage** | `DOMAIN_ONLY` — correct | should never have been recorded at all |

The top-right cell is the one that matters, and it is the worst error this project can make: a
false negative for discovery.

**Therefore the E-value gates whether a hit exists at all, and coverage decides the class among
the hits that survive.** It is not a second classification axis. The class drives one binary
decision — does this reach the bench? — and a `FUNCTIONAL_WEAK` class would either behave like
`FUNCTIONAL`, making it decoration, or like `DOMAIN_ONLY`, in which case it belongs there.

Two consequences:

- **Thresholds are declared per tier, not globally.** In particular there is **no blanket E-value
  floor on T1**: for some short families the gathering threshold is looser than any global cut,
  and overriding curated thresholds would degrade the highest-quality signal in the cascade.
- **`annot_evalue` is the strongest support, not an arbitrary hit's.** Since the best hit is now
  selected by E-value, `annot_evalue` is already the minimum across contributing informative
  hits, so a separate `best_evalue` column would be redundant. `n_informative_hits` is reported
  instead, because 0.9 explained by one domain and 0.9 explained by six fragments are different
  claims.

---

## 7. The artefact screen uses AntiFam's gathering thresholds, not an E-value

Section 1 argues that a curated per-family bit-score cutoff beats a global E-value for
Pfam. The same argument applies with more force to AntiFam, and the size of the effect was
measured rather than assumed.

### The measurement

Every one of AntiFam's **278 profiles carries a GA line**. Each profile also carries its own
`STATS LOCAL FORWARD tau lambda`, which is what HMMER uses to turn a bit score into an
E-value. Converting each curated GA to the E-value it corresponds to at our pinned
`-Z 3,497,616`:

| | |
|---|---|
| profiles with a curated GA | 278 / 278 |
| profiles whose GA is **looser** than E=1e-5 | **274 (98.6%)** |
| median E-value at the curated GA | **7.3e-4** — 73x looser than 1e-5 |
| loosest (`Spurious_ORF_67`, GA 19.2 bits) | E = 0.195 |

Reproduce with:

```python
import math
Z, cur, rows = 3497616, {}, []
for line in open("data/refs/antifam/AntiFam.hmm"):
    if line.startswith("NAME"):            cur = {"name": line.split()[1]}
    elif line.startswith("GA "):           cur["ga"] = float(line.split()[1])
    elif line.startswith("STATS LOCAL FORWARD"):
        _, _, _, tau, lam = line.split();  cur["tau"], cur["lam"] = float(tau), float(lam)
    elif line.startswith("//") and "ga" in cur:
        rows.append(Z * math.exp(-cur["lam"] * (cur["ga"] - cur["tau"]))); cur = {}
print(sum(e > 1e-5 for e in rows), "/", len(rows))
```

### Why it matters more here than anywhere else in the cascade

Everywhere else in the pipeline a missed hit costs an annotation. Here it costs the opposite:
AntiFam exists to catch sequences that are **not proteins** — shadow ORFs on the reverse
complement of a real gene, translated rRNA, repeat-derived ORFs. Those are exactly the things
that pass the whole cascade cleanly, because our selection criterion is "nothing named it" and
nothing names a non-protein. They are also *conserved*, because the real feature underneath
them is conserved, so they survive the multi-lineage test at S9 too.

A missed AntiFam hit is therefore a non-protein sent to the bench with a clean record.
Overriding the curator on 98.6% of the database in that stage is not a defensible default.

### The regression test

`tests/test_scripts_smoke.py::test_the_artefact_screen_uses_antifams_curated_thresholds`
runs the real screen against the real AntiFam over a sequence built from the
`Spurious_ORF_67` consensus and degraded until it sits between the two thresholds: domain
score 32.4 bits, which clears that family's GA of 19.2, and i-Evalue 1.5e-05, which does not
clear a 1e-5 floor. Under `--cut_ga` it is flagged; under `-E 1e-5` it is not.

`antifam_max_evalue` is now `null`, which travels through the same
`cascade.passes_significance` gate the tiers use and means "the curated threshold decided
this". The two gates cannot drift apart because they are the same function.
