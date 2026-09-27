"""eggNOG-mapper over the proteins the cascade named.

Only proteins classed FUNCTIONAL or DOMAIN_ONLY are annotated: the terms describe the
neighbours of dark ORFs for the context features (plasmidann.orthology). Proteins
PlasmidScope annotated take PlasmidScope's own eggNOG-mapper result; the rest are sent to
eggNOG-mapper. Every named protein gets a row; one eggNOG-mapper could not place gets an
empty row.
"""
import csv
import pathlib
import subprocess

import _ctx  # noqa: F401

from plasmidann.orthology import parse_annotations

cfg = snakemake.params.orthology
outdir = pathlib.Path(snakemake.output[0]).parent / "emapper"
outdir.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------------------------
# The query set: proteins the cascade named. A dark protein has nothing for eggNOG to
# transfer an ortholog from, and searching 3.5M sequences to learn that would cost days.
# ------------------------------------------------------------------------------------
named = set()
with open(snakemake.input.prot, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("functional_class") in ("FUNCTIONAL", "DOMAIN_ONLY"):
            named.add(r["seq_id"])

# Proteins PlasmidScope annotated already carry eggNOG-mapper's result from the PlasmidScope
# import; only the other named proteins are sent to eggNOG-mapper here.
from_ps = {}
with open(snakemake.input.ps, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["ps_class"] == "ANNOTATED" and r["seq_id"] in named:
            from_ps[r["seq_id"]] = r

query = outdir / "named.faa"
n_query = 0
with open(query, "w") as out:
    emit = False
    for line in open(snakemake.input.faa):
        if line[0] == ">":
            sid = line[1:].split()[0]
            emit = sid in named and sid not in from_ps
            n_query += emit
        if emit:
            out.write(line)

# When the stage is required, a missing database fails pre-flight. When it is not, the
# search is skipped, every term is recorded as absent, and the columns still exist for the
# context features.
required = bool(cfg.get("required"))
data_dir = pathlib.Path(cfg["data_dir"])
records = {}
if n_query and not data_dir.is_dir():
    if required:
        raise SystemExit(
            f"orthology.required is true but the eggNOG data directory is missing: "
            f"{data_dir}. Run `download_eggnog_data.py --data_dir {data_dir}`, or set "
            "orthology.required to false in config/targets.yaml to run without KEGG terms.")
    print(f"orthology: eggNOG data directory absent ({data_dir}) and not required - "
          "every term recorded as absent, NOT as searched and not found")
elif n_query:
    # DIAMOND mode rather than HMMER mode: this is a bulk ortholog transfer over millions of
    # sequences, and the published defaults are the citable settings (Cantalapiedra et al.
    # 2021, Mol Biol Evol 38:5825). check=True because an empty annotation file here would
    # silently remove the KEGG axis from every context feature downstream.
    subprocess.run(
        f"emapper.py -i {query} -o named --output_dir {outdir} -m diamond "
        f"--cpu {snakemake.threads} --data_dir {cfg['data_dir']} --override",
        shell=True, check=True, stdout=subprocess.DEVNULL)
    annotations = outdir / "named.emapper.annotations"
    if not annotations.exists():
        raise SystemExit(
            f"emapper.py exited 0 but wrote no annotations to {annotations}. The KEGG axis "
            "of every context feature would be silently empty.")
    records = parse_annotations(annotations.read_text())

cols = ["seq_id", "cog_category", "kegg_pathways", "preferred_name", "eggnog_description",
        "eggnog_ogs", "pfams", "gos", "ec", "kegg_ko", "orthology_source"]
with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for sid in sorted(named):
        if sid in from_ps:
            # PlasmidScope publishes no preferred name or description, so those stay
            # empty; its COG/OG identifier stands in for eggnog_ogs.
            p = from_ps[sid]
            w.writerow({"seq_id": sid, "cog_category": p["cog_category"],
                        "kegg_pathways": p["kegg_pathways"], "preferred_name": "",
                        "eggnog_description": "", "eggnog_ogs": p["cog_id"],
                        "pfams": p["pfams"], "gos": p["gos"], "ec": p["ec"],
                        "kegg_ko": p["kegg_ko"], "orthology_source": "plasmidscope"})
            continue
        r = records.get(sid, {})
        w.writerow({"seq_id": sid,
                    "cog_category": r.get("cog_category", ""),
                    "kegg_pathways": ",".join(r.get("kegg_pathways", [])),
                    "preferred_name": r.get("preferred_name", ""),
                    "eggnog_description": r.get("description", ""),
                    "eggnog_ogs": r.get("eggnog_ogs", ""),
                    "pfams": ",".join(r.get("pfams", [])),
                    "gos": ",".join(r.get("gos", [])),
                    "ec": ",".join(r.get("ec", [])),
                    "kegg_ko": ",".join(r.get("kegg_ko", [])),
                    "orthology_source": "emapper" if r else ""})

n_kegg = sum(1 for r in records.values() if r["kegg_pathways"])
n_symbol = sum(1 for r in records.values() if r["preferred_name"])
print(f"orthology: from_plasmidscope={len(from_ps)} queried={n_query} "
      f"annotated={len(records)} with_kegg={n_kegg} with_gene_symbol={n_symbol}")
