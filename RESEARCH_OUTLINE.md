# Research Outline — Environmental spread and architecture of plasmids in IMG/PR
### Extending de Quinto et al. (2026) on multireplicon and small plasmids

*Drafted 2026-06-17. References in `papers_of_interest.txt`. Citations marked `[?]` are from memory — confirm DOIs before manuscript use.*

---

## 0. TL;DR — the one-sentence thesis

> de Quinto et al. derived "rules" for multireplicon plasmid formation from **~24,000 non-redundant, predominantly cultured/clinical** plasmids. We use **IMG/PR's ~700,000 plasmids — including ~30× more, and crucially the *uncultured/environmental* fraction (metagenomes, metatranscriptomes, MAGs, SAGs) with geographic + ecosystem metadata** — to (i) stress-test those rules outside the clinic, (ii) draw the first environment-resolved **biogeography of plasmid replicon architectures**, and (iii) characterise the under-sampled world of **small plasmids**.

The novelty is not "more data." It is that IMG/PR can answer three questions an isolate-derived dataset **structurally cannot**: *where* plasmids live, whether architecture rules are **clinical artefacts vs. universal**, and what the **small/uncultured majority** looks like.

---

## 1. What de Quinto et al. established (the baseline we improve on)

From the preprint (abstract + design):

| Finding | Detail |
|---|---|
| Prevalence | >30% of 24,000 non-redundant plasmids encode **multiple replicons** |
| Phenotype of multireplicon plasmids | Larger; enriched in **AMR + metal + biocide resistance + virulence**; higher **mobility**; broader **host range** |
| Assembly is non-random | Specific replicon **pairs repeatedly fuse**; others rarely fuse despite co-residence |
| Architecture | Replicon pairs sit either in **close proximity** or at **opposite poles** of the plasmid |
| Two evolutionary classes | **Long-term coevolving** pairs vs **transient** associations |
| Mechanism | **Insertion sequences (IS)** drive formation and maintenance |

**Built-in limitations (= our openings):**
1. **Culture & clinical bias.** A non-redundant curated set is dominated by cultivable, clinically/agriculturally relevant taxa. AMR/virulence enrichment may be partly a *sampling* signal.
2. **Reference-based replicon typing is blind to the novel.** Tools like PlasmidFinder `[?]` are trained on Enterobacteriaceae-centric reference replicons → environmental replicons are systematically **untypeable**.
3. **No spatial/ecological axis.** Isolate metadata rarely carries usable coordinates/ecosystem → "where do these plasmids spread?" is unanswerable in their design.
4. **Small plasmids under-represented.** Short, low-copy or non-replicon-typeable plasmids are frequently lost in isolate genome assemblies and from curated DBs.

## 2. What IMG/PR adds (the lever)

699,973 plasmids (154,680 complete) from isolate genomes **and** 28,865 metagenomes, 7,258 metatranscriptomes, 10,499 MAGs, 4,342 SAGs. Per-plasmid metadata we will lean on:

- **Geography**: coordinates → biogeography.
- **Ecosystem**: host-associated / environmental (soil, marine, freshwater, engineered) → ecological niche.
- **Host taxonomy** → host-range breadth.
- **Mobility**: relaxase, T4SS, *oriT* → conjugative vs mobilizable vs non-mobilizable.
- **Resistance**: ARG presence/mechanism.
- **Function**: Pfam / COG / TIGRFAM / KEGG.
- **Completeness & topology**: DTR/ITR → which sequences are safe for architecture analysis.
- **PTUs** (214,950 Plasmid Taxonomic Units; AF≥50%, ANI≥70%) → dereplication unit and host-range proxy.
- **Metatranscriptome-derived plasmids** → *in situ* transcriptional activity (unique to IMG/PR).

---

## 3. Research questions

Each RQ: **Hypothesis → Why it matters → Approach (IMG/PR fields + method) → Relation to de Quinto.**

### Theme A — Do the multireplicon "rules" hold beyond the cultured/clinical world?

**A1. Is the >30% multireplicon prevalence reproduced in metagenome-derived plasmids?**
- *Hypothesis:* Prevalence is **lower** in environmental plasmids (less selection for large composite AMR backbones) but non-trivial.
- *Why:* If >30% is largely clinical, it is a selection signature, not a universal rule of plasmid biology.
- *Approach:* Restrict to **complete** plasmids (154,680). Type replicons with PlasmidFinder `[?]` for known reps **plus** a de-novo rep-protein HMM/clustering layer to capture novel replicons; **report the untypeable fraction explicitly**. Compute multireplicon prevalence stratified by source (isolate / MAG / metagenome) and ecosystem. Dereplicate by PTU before computing prevalence to remove clone over-counting.
- *vs de Quinto:* Direct reproduction + stratification they could not do.

**A2. Are the same replicon pairs the dominant fusions in the environment?**
- *Hypothesis:* A core set of fusions is universal; a second set is **environment-specific**.
- *Why:* Distinguishes intrinsic replicon compatibility from niche-driven assembly.
- *Approach:* Replicon co-occurrence/fusion network per ecosystem; compare edge sets to de Quinto's "repeatedly-fusing" pairs (Jaccard / network alignment). Test their "coexist-but-rarely-fuse" exclusions across ecosystems.
- *vs de Quinto:* Tests universality of their non-random assembly rule.

**A3. Is the AMR/metal/biocide/virulence + host-range enrichment a clinical artefact?**
- *Hypothesis:* Metal/biocide resistance enrichment **persists** environmentally; antibiotic-resistance enrichment **weakens** outside host-associated ecosystems.
- *Why:* Separates anthropogenic-selection cargo from intrinsic multireplicon biology.
- *Approach:* Compare ARG / metal / biocide / virulence / defense gene density (multi vs single replicon) **within each ecosystem** (avoids Simpson's paradox). Host-range = taxonomic breadth of the PTU.
- *vs de Quinto:* Decomposes their pooled enrichment by ecology.

**A4. Is IS-driven formation universal?**
- *Hypothesis:* IS-transposase enrichment at replicon junctions holds across ecosystems but with **different IS families** dominating per biome.
- *Approach:* Pfam transposase density flanking replicons in complete multireplicon plasmids, by ecosystem.
- *vs de Quinto:* Tests generality of their mechanistic claim.

### Theme B — Environmental spread & biogeography (the dimension de Quinto's data cannot reach)

**B1. Replicon-type × ecosystem map — generalists vs specialists.**
- *Hypothesis:* A minority of replicon types are ecological generalists (span many ecosystems); most are specialists.
- *Why:* Generalist replicons are the backbone of cross-ecosystem gene flow ("HGT highways").
- *Approach:* Replicon × ecosystem occurrence matrix; niche-breadth index (e.g. Levins' B / Shannon over ecosystems) per replicon type. **First figure of the project.**
- *vs de Quinto:* Net-new axis.

**B2. Are multireplicon plasmids ecological generalists?**
- *Hypothesis:* Multireplicon plasmids occupy **broader ecosystem ranges**, mirroring their broader host range.
- *Why:* Tests whether "second replicon = expansion device" extends from host range to *habitat* range.
- *Approach:* Ecosystem-breadth (multi vs single replicon), controlling for plasmid size and sampling effort.

**B3. Biogeography of PTUs.**
- *Hypothesis:* Some PTUs show biome/latitudinal structuring; others are cosmopolitan.
- *Why:* Cosmopolitan PTUs imply global dissemination routes; structured PTUs imply local evolution.
- *Approach:* PTU occurrence vs coordinates; spatial autocorrelation (Moran's I) and biome enrichment; rarefy by sampling effort per region.

**B4. Environmental reservoirs of clinically relevant / AMR plasmids (One Health).**
- *Hypothesis:* Clinically dominant replicon types and AMR-carrying multireplicon plasmids have identifiable **environmental reservoirs** (e.g. wastewater, agricultural soil, freshwater).
- *Why:* The highest-impact human-health extension of de Quinto's AMR focus — *where* resistance backbones live before/after the clinic.
- *Approach:* For clinical replicon types and ARG-positive plasmids, map non-clinical ecosystem occurrence; flag PTUs bridging clinical and environmental samples.
- *vs de Quinto:* Turns their AMR-spread finding into a spatial/reservoir map.

### Theme C — Small plasmids: the under-explored majority

**C1. Size spectrum and the small-plasmid mode.**
- *Hypothesis:* IMG/PR's metagenomic capture reveals a far larger small-plasmid (<10 kb) fraction than isolate-only DBs; many are single-replicon or rolling-circle (replicon-typeable only by Rep_2/Rep_trans HMMs).
- *Approach:* Length distribution by source; quantify small/complete plasmids; classify mobility (mobilizable vs non-mobilizable vs cryptic).

**C2. What do small plasmids carry — cryptic, or functional?**
- *Hypothesis:* A substantial small-plasmid fraction carries **defense systems (anti-phage), toxin–antitoxin, restriction–modification, or single resistance genes**, not "nothing."
- *Why:* Reframes small plasmids from molecular curiosities to mobile functional modules; the anti-phage-defense angle `[?]` (Rocha & Bikard 2022) is timely and largely unexplored at IMG/PR scale.
- *Approach:* Functional density (Pfam/KO) of small vs large plasmids; defense-system annotation (e.g. PADLOC/DefenseFinder-style HMMs); fraction with no detectable cargo ("truly cryptic").

**C3. Do small plasmids hitchhike?**
- *Hypothesis:* Mobilizable small plasmids co-occur in samples/hosts with conjugative (often multireplicon) plasmids that supply transfer machinery.
- *Approach:* Same-sample/same-host co-occurrence of mobilizable-non-conjugative small plasmids with conjugative plasmids; relaxase-vs-T4SS complementation logic.

**C4. Are small-plasmid replicons the building blocks of multireplicon backbones? (links C→A)**
- *Hypothesis:* Replicons that appear standalone in small plasmids recur as **components** of multireplicon plasmids.
- *Why:* A mechanistic bridge: small plasmids as the monomers that fuse into composite backbones.
- *Approach:* Cross-reference replicon types between the small-plasmid set and the multireplicon component set; test enrichment.

### Theme D — New perspectives (higher-risk, higher-novelty)

**D1. Plasmid "guild" ecology.** Bipartite network (replicons × ecosystems); community-detect ecological guilds of plasmids — an ecological complement to taxonomic PTUs.

**D2. Anti-phage defense as an alternative spread driver.** Test whether **defense-system cargo** predicts host/ecosystem breadth as well as, or better than, AMR cargo. Reframes the "what makes a plasmid spread" narrative beyond resistance. `[?]`

**D3. Second replicon as a host-range expander (directional test).** Within PTUs, test whether multireplicon members reach broader host taxonomy than single-replicon relatives — evidence that replicon acquisition *accompanies* host jumps.

**D4. Fusion intermediates at population scale.** Use completeness + IS-junction signatures to identify putative **fusion/deletion intermediates** — population-scale snapshots of the fusion-deletion life cycle described in the companion 2026 bioRxiv plasmid-life-cycle preprint.

**D5. In-situ activity from metatranscriptomes.** Which replicon types / cargo classes are **transcriptionally active** in metatranscriptome-derived plasmids? Impossible in any isolate DB; uniquely enabled by IMG/PR.

---

## 4. Methodological caveats (read before any claim)

1. **Completeness/assembly bias** — restrict architecture (multireplicon, IS-junction, length) analyses to **complete** plasmids; use fragments only for coarse prevalence, flagged.
2. **Replicon-typing blind spots** — PlasmidFinder is reference-biased; supplement with Rep-protein HMMs / de-novo clustering and **always report the untypeable %**. "No detected replicon" ≠ "no replicon."
3. **Redundancy & sampling effort** — dereplicate to PTU (or ANI) before prevalence/enrichment stats; weight/rarefy by per-ecosystem and per-region sampling effort. Clinical clones are massively over-sampled.
4. **Metadata completeness** — quantify how many plasmids actually have coordinates/ecosystem/host **before** biogeography claims; report coverage.
5. **Host prediction uncertainty** — metagenomic host assignment is probabilistic; propagate confidence, don't treat as ground truth.
6. **Correlation ≠ causation** — "second replicon ↔ broader range" is associational; frame as hypothesis-generating, note experimental follow-up (the de Quinto group's wet-lab strength).
7. **Cross-DB comparability** — to compare directly to de Quinto's 24k, re-type their plasmids and IMG/PR with the *same* pipeline; don't compare across different replicon-typing schemes.

## 5. First milestones (verifiable)

1. **Download IMG/PR** → `data/img_pr/` → verify ≈699,973 sequences + metadata rows. *Verify:* `seqkit stats` count matches metadata row count.
2. **Reproduce baseline** → multireplicon prevalence on IMG/PR complete plasmids. *Verify:* falls in a defensible range vs de Quinto's 30% (and explainable if it differs).
3. **Theme B1 figure** → replicon × ecosystem matrix + niche-breadth ranking. *Verify:* generalist replicons are biologically sensible (e.g. broad-host-range Inc groups rank high).
4. **Theme C1 census** → small-plasmid size/mobility/function table. *Verify:* counts reconcile with total; size modes are reproducible across sources.

## 6. Open decisions for you

- **Replicon-typing strategy** — PlasmidFinder-only (comparable to de Quinto, but misses environmental reps) vs PlasmidFinder + de-novo HMM layer (more complete, more work). *Recommendation: both, reporting untypeable fraction.*
- **Dereplication unit** — use IMG/PR's PTUs as-is, or recompute clustering for full control. *Recommendation: start with IMG/PR PTUs; recompute only if needed.*
- **Scope of first paper** — narrow (Theme B environmental spread, the cleanest novel story) vs broad (A+B+C). *Recommendation: lead with B + C1/C2; keep A as the validation/baseline; D as follow-ups.*
