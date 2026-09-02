#!/usr/bin/env python3
"""LOCKED environment reconciliation (2026-07-10) — the canonical producer of the habitat taxonomy.

Maps the raw per-plasmid environment labels (IMG/PR GOLD + NCBI BioSample + mMGE SRA) onto one locked
two-tier taxonomy, plus three orthogonal tags. Supersedes the earlier notebook-inline "PROPOSED" version.

Locked decisions (2026-07-10):
  * Environment (habitat) is THE label. Top level = 3 GOLD-aligned domains: Host-associated /
    Environmental / Engineered. Non-habitat buckets: `Simulated-artifact` (GOLD Modeled;Simulated
    communities, excluded from analysis) and `Unknown` (unresolved + unlabelled + geography-only merged).
  * Human habitats keep body-site granularity (blood / urine / respiratory / skin-wound / oral /
    gut-faeces / other-clinical / unspecified) — this is how we know "human feces" vs "wastewater" etc.
  * `is_clinical` = additional tag (NOT a top level): clear clinical body site OR PLSDB disease tag OR
    hospital/nosocomial/patient context. Environment-label-based, NOT source-DB based.
  * `lifestyle` / `sample_nature` (the metagenomic-vs-isolate axis) are DROPPED — the source-DB tag does
    not reflect true sample origin (clinical datasets contain metagenomic samples). Database provenance
    lives in the `sources` column of the working set instead.
  * Ambiguous bare `feces/stool` -> Human:gut/faeces; bare `gut` -> unspecified host.

Inputs:  data/plasmidscope_primary/{analysis_table.tsv, plsdb_enrichment.tsv}
Output:  data/plasmidscope_primary/environment_reconciled.tsv
Run:     python3 scripts/reconcile_environment.py
"""
import csv
import re
import sys
from collections import Counter

csv.field_size_limit(sys.maxsize)
PP = "data/plasmidscope_primary"
ANALYSIS = f"{PP}/analysis_table.tsv"
PLSDB = f"{PP}/plsdb_enrichment.tsv"
OUT = f"{PP}/environment_reconciled.tsv"

R = lambda p: re.compile(p, re.I)
# ordered free-text rules: FIRST match wins (specific/high-priority first; wastewater before water)
FREETEXT_RULES = [
    (R(r"waste ?water|sewage|activated sludge|\bsludge\b|wwtp|influent|effluent|treatment plant"), "Engineered", "Wastewater/sewage"),
    (R(r"bioreactor|fermentor|fermenter|anaerobic digest"), "Engineered", "Bioreactor"),
    (R(r"kimchi|dairy|fermentation|\bcheese\b|\byogurt|\bmilk\b|\bmeat\b|\bpork\b|\bfood\b|sausage|sauerkraut"), "Engineered", "Food/fermentation"),
    (R(r"built environment|hospital surface|indoor|building|door handle|\bsink\b|catheter|\bplastic\b"), "Engineered", "Built environment"),
    (R(r"compost|solid waste|landfill|manure|biosolid"), "Engineered", "Solid waste/compost"),
    (R(r"bioremediat|contaminated|oil spill|petroleum|hydrocarbon|industrial"), "Engineered", "Industrial/remediation"),
    (R(r"blood|bacteremia|bacteraemia|bloodstream|\bserum\b|\bplasma\b"), "Host-associated", "Human: blood"),
    (R(r"urine|urinary|\buti\b|bladder|catheter urine"), "Host-associated", "Human: urine"),
    (R(r"sputum|respiratory|\blung|nasal|nasopharyn|bronch|throat|\bcough|pneumon|trachea"), "Host-associated", "Human: respiratory"),
    (R(r"\bskin\b|\bwound|\bpus\b|abscess|ulcer|\bsurgical site"), "Host-associated", "Human: skin/wound"),
    (R(r"\boral\b|saliva|\bdental|\btooth|\bgum\b|periodont|oral cavity|mouth"), "Host-associated", "Human: oral"),
    (R(r"cerebrospinal|\bcsf\b|\bbile\b|periton|\bear\b|\beye\b|conjunctiv|vagina|cervic|genital|\bsemen"), "Host-associated", "Human: other clinical"),
    (R(r"human (gut|stool|fecal|faecal|feces|faeces|intest)|human gut|gut metagenome|gut microbiome|rectal|\bstool\b|\bfeces\b|\bfaeces\b|colon|intestin|\bileu|gastro"), "Host-associated", "Human: gut/faeces"),
    (R(r"\bhuman\b|homo sapiens|patient|clinical|hospital|nursing home|healthy (person|volunteer|college|adult)|nosocomial"), "Host-associated", "Human: unspecified"),
    (R(r"\bpig\b|\bswine\b|\bsus scrofa|porcine|piglet|\bboar\b"), "Host-associated", "Livestock: pig"),
    (R(r"\bcattle\b|\bcow\b|\bbovine\b|\bbos taurus|\bcalf\b|\bbeef\b|\bveal\b|dairy cow"), "Host-associated", "Livestock: cattle"),
    (R(r"\bsheep\b|\bovine\b|\bgoat\b|\blamb\b|caprine|\bcamel\b|\bhorse\b|equine|\bmule\b"), "Host-associated", "Livestock: other"),
    (R(r"chicken|\bpoultry\b|\bgallus\b|\bturkey\b|\bduck\b|\bgoose\b|\bbroiler\b|\bhen\b|\bavian\b|\bbird\b"), "Host-associated", "Poultry/bird"),
    (R(r"\bdog\b|\bcat\b|canis|feline|felis|companion animal|\bpet\b"), "Host-associated", "Companion animal"),
    (R(r"\brat\b|rattus|\bmouse\b|\bmice\b|\bmurine\b|rodent|\bvole\b|\bhamster\b"), "Host-associated", "Rodent"),
    (R(r"salmo|oncorhynchus|\bfish\b|\btrout\b|\bsalmon\b|penaeus|\bshrimp\b|\bprawn\b|aquaculture|\btilapia\b|\bcarp\b"), "Host-associated", "Fish/aquaculture"),
    (R(r"ixodes|\btick\b|nasonia|\binsect|drosophila|\bbee\b|\baphid\b|\bmosquito\b|\bfly\b|arthropod|\blouse\b|\bmoth\b"), "Host-associated", "Insect/arthropod"),
    (R(r"root nodule|rhizosphere|phyllosphere|\bplant\b|\bleaf\b|\broot\b|\bseed\b|phaseolus|trifolium|medicago|glycine|oryza|\btomato\b|\bmaize\b|\bwheat\b|\brice\b|malus|\bwalnut\b|legume|\bgrass\b|\bcrop\b"), "Host-associated", "Plant"),
    (R(r"nodule|\bleaves\b|\bleaf\b|phyllosph|infected tissue|\bcitrus\b|\bpotato\b|floral|\bfruit\b|\bflower\b|\bstem\b|\bbark\b"), "Host-associated", "Plant"),
    (R(r"skeletonema|\balgae\b|\bdiatom\b|phytoplankton|cyanobacter"), "Host-associated", "Algae"),
    (R(r"\bfungi\b|\bfungal\b|\byeast\b|mycel"), "Host-associated", "Fungi"),
    (R(r"porifera|\bsponge\b|mollusc|\bcoral\b|annelid|\bworm\b|invertebrat|nematod|cephalochord"), "Host-associated", "Other invertebrate"),
    (R(r"mammal|\bmilk\b"), "Host-associated", "Other mammal"),
    (R(r"\bgut\b|\bfecal\b|\bfaecal\b|\bfeces\b|\bfaeces\b|cloacal|intestin|\bcaecal\b|\bcecal\b|\brumen\b|manure"), "Host-associated", "Gut/faeces: unspecified host"),
    (R(r"hydrothermal|\bvent\b|hot spring|geothermal|\bsalt\b|\bbrine\b|hypersaline|\bsaline\b|sea ice|\bglacier|permafrost|\bsnow\b"), "Environmental", "Extreme/other"),
    (R(r"agricultural|farmland|\bfield\b|\bcrop field|pasture|\bmeadow\b|\bforest\b|\bgrassland\b"), "Environmental", "Terrestrial/soil"),
    (R(r"seawater|sea water|\bmarine\b|\bocean\b|estuar|\bcoastal\b|\bbrackish\b"), "Environmental", "Aquatic: marine"),
    (R(r"freshwater|river|\blake\b|\bpond\b|\bstream\b|\bdrinking water|groundwater|\bwater\b|sediment|\baquatic\b|\bwetland\b|\blagoon\b"), "Environmental", "Aquatic: freshwater/other"),
    (R(r"\bsoil\b|\bterrestr|\bsand\b|\bmud\b|\bdust\b|\bcave\b|\brock\b"), "Environmental", "Terrestrial/soil"),
    (R(r"\bair\b|\baerosol\b|atmospher"), "Environmental", "Air"),
    (R(r"\benvironment\b|\benvironmental\b|\bmetagenome\b|\bbiofilm\b"), "Environmental", "Environmental: unspecified"),
]
CLINICAL_SITES = {"Human: blood", "Human: urine", "Human: respiratory", "Human: skin/wound", "Human: other clinical"}
CLINICAL_RX = R(r"\bclinical\b|\bpatient\b|hospital|nosocomial|bacter[ae]mia|\bicu\b|\bsepsis\b|\bwound\b|blood culture")
TOP_ORDER = ["Host-associated", "Environmental", "Engineered", "ARTIFACT"]


def classify_freetext(text):
    for rx, top, sub in FREETEXT_RULES:
        if rx.search(text):
            return top, sub
    return None


def classify_gold(path):
    parts = [p.strip() for p in path.split(";")]
    l1, l2 = parts[0], parts[1] if len(parts) > 1 else ""
    if l1.startswith("Engineered") and l2 == "Modeled":
        return "ARTIFACT", "Simulated communities"
    if l1.startswith("Engineered"):
        # lab enrichment cultures + content-free bare 'Engineered' are NOT a real environmental-
        # engineering site (no field habitat) -> excluded like the simulated-community artifact.
        if l2 == "Lab enrichment":
            return "LAB-ARTIFACT", "Lab enrichment"
        if not l2:
            return "LAB-ARTIFACT", "Engineered: other"
        return "Engineered", {"Wastewater": "Wastewater/sewage", "WWTP": "Wastewater/sewage",
            "Bioreactor": "Bioreactor", "Built environment": "Built environment", "Solid waste": "Solid waste/compost",
            "Food production": "Food/fermentation", "Bioremediation": "Industrial/remediation",
            "Industrial production": "Industrial/remediation", "Biotransformation": "Industrial/remediation",
            "Artificial ecosystem": "Built environment"}.get(l2, l2)
    if l1.startswith("Environmental"):
        return "Environmental", {"Aquatic": "Aquatic: freshwater/other", "Terrestrial": "Terrestrial/soil", "Air": "Air"}.get(l2, l2 or "Environmental: unspecified")
    if l1.startswith("Host-associated"):
        return "Host-associated", {"Mammals: Human": "Human: unspecified", "Mammals": "Other mammal",
            "Arthropoda: Insects": "Insect/arthropod", "Arthropoda: Crustaceans": "Fish/aquaculture", "Plants": "Plant",
            "Algae": "Algae", "Fungi": "Fungi", "Birds": "Poultry/bird", "Fish": "Fish/aquaculture",
            "Annelida": "Other invertebrate", "Porifera": "Other invertebrate", "Mollusca": "Other invertebrate",
            "Invertebrates": "Other invertebrate", "Cephalochordata": "Other invertebrate",
            "Microbial": "Environmental: unspecified"}.get(l2, l2 or "Host: other")
    return None


def reconcile(gold, host, iso, sra, geo):
    if gold:
        cands = [c for c in (classify_gold(p.strip()) for p in gold.split(" | ") if p.strip()) if c]
        real = sorted([c for c in cands if c[0] not in ("ARTIFACT", "LAB-ARTIFACT")],
                      key=lambda c: TOP_ORDER.index(c[0]))
        if real:
            return (*real[0], "GOLD")
        if cands:
            return (*cands[0], "GOLD")
    bio = f"{host} ; {iso}".strip(" ;")
    if bio:
        c = classify_freetext(bio)
        if c:
            return (*c, "BioSample")
    if sra.strip():
        c = classify_freetext(sra)
        if c:
            return (*c, "SRA")
    if geo.strip():
        return ("LOCATION-ONLY", "Geography only (no habitat)", "BioSample")
    return ("UNRESOLVED", "Unresolved", "")


# LOCKED top-level renaming: 3 domains kept; artifact + non-habitat buckets collapsed
LOCK_TOP = {"ARTIFACT": "Simulated-artifact", "LAB-ARTIFACT": "Lab-artifact",
            "LOCATION-ONLY": "Unknown", "UNRESOLVED": "Unknown", "UNLABELLED": "Unknown"}


def main():
    disease = {}
    eco = {}
    for r in csv.DictReader(open(PLSDB, newline=""), delimiter="\t"):
        disease[r["plasmid_id"]] = r.get("plsdb_disease_tags", "").strip()
        eco[r["plasmid_id"]] = r.get("plsdb_ecosystem_tags", "").strip().lower()

    cols = ["plasmid_id", "hab_top", "hab_sub", "hab_channel", "is_clinical"]
    top_c, clin = Counter(), 0
    with open(ANALYSIS, newline="") as fh, open(OUT, "w", newline="") as out:
        r = csv.DictReader(fh, delimiter="\t")
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for row in r:
            rep = row["plasmid_id"]
            if row["env_status"] != "labelled":
                top, sub, chan = "UNLABELLED", "Unlabelled", ""
            else:
                top, sub, chan = reconcile(row["gold_ecosystem"], row["biosample_host"],
                                           row["biosample_isolation_source"], row["sra_sample"],
                                           row["biosample_geo"])
            # clinical tag (additional): clear clinical body site OR PLSDB disease OR clinical free-text/hospital
            blob = f'{row["biosample_isolation_source"]} ; {row["sra_sample"]} ; {row["gold_ecosystem"]}'.lower()
            is_clin = (sub in CLINICAL_SITES or bool(disease.get(rep)) or "hospital" in eco.get(rep, "")
                       or bool(CLINICAL_RX.search(blob)))
            hab_top = LOCK_TOP.get(top, top)
            if is_clin:
                clin += 1
            top_c[hab_top] += 1
            w.writerow({"plasmid_id": rep, "hab_top": hab_top, "hab_sub": sub,
                        "hab_channel": chan, "is_clinical": "1" if is_clin else "0"})

    print(f"wrote {OUT}")
    for k, v in top_c.most_common():
        print(f"  {k:<20} {v:>8,}")
    print(f"  is_clinical=1        {clin:>8,} ({100*clin/sum(top_c.values()):.1f}%)")


if __name__ == "__main__":
    main()
