#!/usr/bin/env python3
"""Phase A.4 + A.6: assemble the PlasmidScope-primary working set.

For every deduplicated complete plasmid, gather ALL raw environment labels across its source channels
(kept verbatim, decision #3/#4): IMG/PR GOLD ecosystem, NCBI BioSample (per INSDC accession), and mMGE
SRA sample. Flag lab-made/synthetic plasmids from the recovered organism (A.4, conservative + logged)
and discard them. Environment is NOT a filter — blanks are kept as `unknown`. Produce the working set
and the final funnel/ledger report.

Inputs (data/plasmidscope_primary/ unless noted):
  complete_provenance.tsv, env_local.tsv, sra_env_raw.tsv
  data/processed/plasmidscope/biosample/env_raw.tsv
Outputs:
  data/plasmidscope_primary/working_set.tsv
  data/plasmidscope_primary/discarded_labmade.tsv
  reports/plasmidscope_primary_working_set.md
Run:  python3 scripts/assemble_working_set.py
"""
import csv
import os
import sys
from collections import Counter

csv.field_size_limit(sys.maxsize)
PP = "data/plasmidscope_primary"
PROV = f"{PP}/complete_provenance.tsv"
ENV_LOCAL = f"{PP}/env_local.tsv"
ENV_RAW = "data/processed/plasmidscope/biosample/env_raw.tsv"
SRA = f"{PP}/sra_env_raw.tsv"
OUT = f"{PP}/working_set.tsv"
DISC = f"{PP}/discarded_labmade.tsv"
REPORT = "reports/plasmidscope_primary_working_set.md"

LABMADE_SIGNALS = ["synthetic construct", "cloning vector", "expression vector", "shuttle vector"]


def strip_ver(a):
    return a.rsplit(".", 1)[0]


def labmade_hit(*texts):
    blob = " ; ".join(t for t in texts if t).lower()
    for s in LABMADE_SIGNALS:
        if s in blob:
            return s
    return ""


def main():
    # BioSample env keyed by version-stripped accession (prefer rows that actually have env)
    bios = {}
    with open(ENV_RAW, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = strip_ver(r["accession"])
            if k not in bios or (r["status"] == "env" and bios[k]["status"] != "env"):
                bios[k] = r
    # SRA env keyed by run
    sra = {}
    if os.path.exists(SRA):
        with open(SRA, newline="") as fh:
            sra = {r["run"]: r for r in csv.DictReader(fh, delimiter="\t")}
    # per-plasmid parsed channels
    local = {}
    with open(ENV_LOCAL, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            local[r["rep_id"]] = r

    cols = ["plasmid_id", "sources", "lifestyle", "topology", "size_bp", "n_raw_entries",
            "gold_ecosystem", "biosample_isolation_source", "biosample_host", "biosample_geo",
            "organism", "sra_sample", "env_status", "disposition"]
    os.makedirs(PP, exist_ok=True)
    kept = discarded = 0
    env_channel = Counter()
    has_any_env = 0
    life_kept = Counter()
    disc_rows = []

    with open(PROV, newline="") as fh, open(OUT, "w", newline="") as out:
        r = csv.DictReader(fh, delimiter="\t")
        w = csv.writer(out, delimiter="\t")
        w.writerow(cols)
        for row in r:
            rep = row["plasmid_id"].split(",")[0]
            lc = local.get(rep, {})
            gold = lc.get("imgpr_ecosystem_raw", "")
            accs = [a for a in lc.get("insdc_accessions", "").split(";") if a]
            runs = [x for x in lc.get("sra_runs", "").split(";") if x]

            # BioSample channel (isolate INSDC accessions): organism here is the bacterial host, NOT
            # a habitat, so it does not count as environment.
            b_iso, b_host, b_geo, b_org = [], [], [], []
            for a in accs:
                b = bios.get(strip_ver(a))
                if b:
                    if b.get("isolation_source"): b_iso.append(b["isolation_source"])
                    if b.get("host"): b_host.append(b["host"])
                    if b.get("geo_loc_name"): b_geo.append(b["geo_loc_name"])
                    if b.get("organism"): b_org.append(b["organism"])
            # SRA channel (mMGE): organism IS the habitat (e.g. 'human gut metagenome'), so it counts.
            s_iso, s_host, s_body, s_org = [], [], [], []
            for run in runs:
                s = sra.get(run)
                if s:
                    if s.get("isolation_source"): s_iso.append(s["isolation_source"])
                    if s.get("host"): s_host.append(s["host"])
                    if s.get("body_site"): s_body.append(s["body_site"])
                    if s.get("organism"): s_org.append(s["organism"])

            uniq = lambda xs: " | ".join(dict.fromkeys(x for x in xs if x))
            gold_s = gold
            biosample_iso, biosample_host, biosample_geo = uniq(b_iso), uniq(b_host), uniq(b_geo)
            sra_sample = uniq(s_iso + s_body + s_org)   # mMGE raw habitat text (comment/body/metagenome)
            org_s = uniq(b_org + s_org)                 # all organisms, for lab-made detection
            meta_org = any("metagenome" in o.lower() for o in s_org)

            biosample_env = bool(biosample_iso or biosample_host or biosample_geo)
            sra_env = bool(s_iso or s_host or s_body or meta_org)
            if gold_s: env_channel["gold"] += 1
            if biosample_env: env_channel["biosample"] += 1
            if sra_env: env_channel["sra"] += 1
            any_env = bool(gold_s or biosample_env or sra_env)
            if any_env: has_any_env += 1
            env_status = "labelled" if any_env else "unknown"

            # A.4 lab-made flag from recovered organism + PlasmidScope host
            sig = labmade_hit(org_s, row.get("host", ""))
            if sig:
                discarded += 1
                disc_rows.append([rep, row["sources"], sig, org_s or row.get("host", "")])
                continue
            kept += 1
            life_kept[row["lifestyle"]] += 1
            w.writerow([rep, row["sources"], row["lifestyle"], row["topology"], row["size_bp"],
                        row["n_raw_entries"], gold_s, biosample_iso, biosample_host, biosample_geo,
                        org_s, sra_sample, env_status, "keep"])

    with open(DISC, "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(["plasmid_id", "sources", "labmade_signal", "organism_or_host"])
        wr.writerows(disc_rows)

    total_complete = kept + discarded
    L, w = [], lambda s: L.append(s)
    w("# PlasmidScope-primary — working set & environment (Phase A.4 + A.6)\n")
    w("**Pipeline:** `scripts/assemble_working_set.py` merging `complete_provenance.tsv`, "
      "`env_local.tsv` (IMG/PR GOLD, `scripts/recover_env_local.py`), the NCBI BioSample table "
      "`env_raw.tsv` (`scripts/fetch_biosample_env.py`), and the mMGE SRA table `sra_env_raw.tsv` "
      "(`scripts/fetch_sra_env.py`). Regenerate:\n```\npython3 scripts/assemble_working_set.py\n```\n")

    w("## Inclusion funnel (full accounting)\n")
    w("| stage | n | note |")
    w("|---|---:|---|")
    w("| PlasmidScope `ALL` (deduplicated) | 852,600 | MMseqs2 100% id/cov |")
    w(f"| complete / closed | {total_complete:,} | −644,240 incomplete/unknown (Phase A.2) |")
    w(f"| − lab-made / synthetic | −{discarded:,} | flagged from organism/host (see below) |")
    w(f"| **= WORKING SET** | **{kept:,}** | everything else kept, incl. `unknown` environment |\n")

    w("## Lab-made discards (conservative, logged)\n")
    w(f"**{discarded:,}** plasmids discarded as lab-made, each logged in `discarded_labmade.tsv` with "
      "the matching signal. Signals used: " + ", ".join(f"`{s}`" for s in LABMADE_SIGNALS) + ".\n")
    if disc_rows:
        sc = Counter(d[2] for d in disc_rows)
        w("| signal | n |")
        w("|---|---:|")
        for s, c in sc.most_common():
            w(f"| {s} | {c:,} |")
        w("")

    w("## Environment coverage of the working set (raw, un-reconciled)\n")
    w(f"- Any environment label recovered: **{has_any_env:,} / {kept:,} "
      f"({100*has_any_env/max(kept,1):.1f}%)**; the rest are `unknown` (kept, decided later).")
    w("| channel | plasmids with a label from it |")
    w("|---|---:|")
    w(f"| IMG/PR GOLD ecosystem | {env_channel['gold']:,} |")
    w(f"| NCBI BioSample (INSDC) | {env_channel['biosample']:,} |")
    w(f"| mMGE SRA sample | {env_channel['sra']:,} |")
    w("\n(A plasmid can draw from several channels — provenance and environments are kept as sets.)\n")

    w("## Working set by lifestyle\n")
    w("| lifestyle | n | % |")
    w("|---|---:|---:|")
    for k in ("metagenomic", "isolate", "mixed", "other"):
        if life_kept.get(k):
            w(f"| {k} | {life_kept[k]:,} | {100*life_kept[k]/kept:.1f}% |")
    w("\n**Environment reconciliation (raw → harmonised habitat taxonomy) is deferred**, per the "
      "locked plan — labels are kept verbatim here for inspection first.\n")

    with open(REPORT, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"wrote {OUT} (kept {kept:,}), {DISC} (discarded {discarded:,}), {REPORT}")
    print(f"any_env={has_any_env:,} ({100*has_any_env/max(kept,1):.1f}%) "
          f"gold={env_channel['gold']:,} biosample={env_channel['biosample']:,} sra={env_channel['sra']:,}")


if __name__ == "__main__":
    main()
