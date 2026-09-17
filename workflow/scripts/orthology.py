"""S4b: eggNOG-mapper over the proteins the cascade named.

Runs on the ANNOTATED fraction, not the dark one. See plasmidann.orthology for why a
dark-protein pipeline spends compute describing the known genes: S8 asks what a dark ORF's
neighbours do, the cascade answers in free text, and free text cannot be aggregated into
pathways. FESNov's neighbourhood metric is defined over KEGG pathway membership because a
pathway is a term you can count.

Nothing is filtered. Every unique protein gets a row; the ones eggNOG could not place get
an empty one, which is the honest record of an absent term rather than a missing row.
"""
import _ctx  # noqa: F401
import csv
import pathlib
import subprocess

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

query = outdir / "named.faa"
n_query = 0
with open(query, "w") as out:
    emit = False
    for line in open(snakemake.input.faa):
        if line[0] == ">":
            emit = line[1:].split()[0] in named
            n_query += emit
        if emit:
            out.write(line)

# Same contract as structural evidence: when the stage is declared required, a missing
# database fails PRE-FLIGHT in seconds rather than here. When it is not required the stage
# is skipped and every term is recorded as absent - which is honest - and the columns still
# exist, because S8 reads them.
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
            "of every S8 context feature would be silently empty.")
    records = parse_annotations(annotations.read_text())

cols = ["seq_id", "cog_category", "kegg_pathways", "preferred_name", "eggnog_description",
        "eggnog_ogs", "pfams", "gos", "ec", "kegg_ko"]
with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for sid in sorted(named):
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
                    "kegg_ko": ",".join(r.get("kegg_ko", []))})

n_kegg = sum(1 for r in records.values() if r["kegg_pathways"])
n_symbol = sum(1 for r in records.values() if r["preferred_name"])
print(f"orthology: queried={n_query} annotated={len(records)} with_kegg={n_kegg} "
      f"with_gene_symbol={n_symbol}")
