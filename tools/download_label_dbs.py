"""Download the protein label databases into data/refs/labels/<db>/.

WHAT IS INSTALLED

  tadb      TADB 3.0, experimentally validated entries (*_exp.fas), protein records only
  bacmet    BacMet 2.0, experimentally confirmed genes (BacMet2_EXP)
  oritdb    oriTDB 2.0, relaxase, auxiliary protein and T4CP (validated plus predicted)
  card      CARD, the current release: card.json and the protein homolog model FASTA
  mobileog  mobileOG-db beatrix-1.6, every entry (Manual, Homology, Keyword Search)
  dbapis    dbAPIS, the verified APIS proteins and their sequence homologues
  acrdb     Anti-CRISPRdb version 2.2, every entry of the core dataset

EVIDENCE RULE

mobileOG-db, dbAPIS and Anti-CRISPRdb keep all evidence classes, and the class of every
entry is recorded in its FASTA header, so each label can carry it.
TADB keeps its experimentally validated files (*_exp) and BacMet its experimentally
confirmed set (BacMet2_EXP); CARD keeps every protein homolog model; oriTDB keeps
validated plus predicted entries. Each SOURCE file states its rule and the counts.

Each directory holds the downloads under raw/, the files built from them, VERSION
(one line), SOURCE (url, citation, licence, download date, FASTA header format) and
MANIFEST.sha256 (every file, checkable with `sha256sum -c`).

PINNED DOWNLOADS

Every download except CARD is pinned by SHA-256, so a changed file upstream stops the
install instead of silently changing the labels. CARD is taken as the current release, and
its version is read from the download (card.json). A pinned file may instead be taken from a local copy (--local DIR, looked up as
DIR/<db>/<file name>), and only when its SHA-256 matches.

IDEMPOTENT

A database whose VERSION, SOURCE and MANIFEST.sha256 are present, whose recorded version
is the expected one and whose files all match the manifest is skipped, so a second run
downloads nothing.

Run from the repository root:  python tools/download_label_dbs.py [--local DIR ...]
"""
import argparse
import csv
import datetime
import hashlib
import json
import pathlib
import shutil
import sys
import tarfile
import urllib.request

TADB_URL = "https://bioinfo-mml.sjtu.edu.cn/TADB3/download/"
TADB_FILES = {  # the experimentally validated files listed on TADB3/download.php
    "regulator_exp.fas": "a52d85add064f4825b4ec14bb97b440f55ec40e688fa1bbf2877bf902a2d8427",
    "type_I_AT_exp.fas": "ca881c9391cf61e5b42d96bf60a97ce84fc5ee1370983ab1057d7e704ed2ee86",
    "type_I_T_exp.fas": "d51e00f7546a7905b0775becd467aa0797793d99de2046d78542b3b502375dd7",
    "type_II_AT_exp.fas": "07dac6cfcfad815a4d2c92f5487e0ad327fa71a2fe558f9ffa319a7a081816c5",
    "type_II_T_exp.fas": "3e9166ffc2fb8a301597c3b4671eb7f610afb24550d7345d1961fd05f59fd865",
    "type_III_AT_exp.fas": "d51f18be2ec889b5c9cabb1a60b901b080b4db9dd9269621c7f420af720d40f8",
    "type_III_T_exp.fas": "4a78b59ba31b395e952d80bb7004b67dbe0f12ba5618e4cc770e1b9bd76b965c",
    "type_IV_AT_exp.fas": "3231cab7740edc7cc01c2e0845c5dee9df05560066788187acb789a08f6ac68c",
    "type_IV_T_exp.fas": "d7b72d43809b058d383e8ae0d4393ac19b9b0c6efb659f4e5dcf8ca8e5b596eb",
    "type_V_AT_exp.fas": "4027fdee8c2b9a309320fd1ea70a1b4148f149ffd8615623b9413bc43814075c",
    "type_V_T_exp.fas": "6d94156ffee408bd2ab6186a57ea0944690b16f99b0a0a369919adea276c9e90",
    "type_VI_AT_exp.fas": "48e1e6e75743bdd8b440fadccbca91b6eea14d7556c1590bbfeb445350ec2420",
    "type_VI_T_exp.fas": "dc667130adb5f700b58788ac548a05163cca14fc34b0980205c738a69b2735a9",
    "type_VII_AT_exp.fas": "a75a90fc1740fb82e3d39166956a3e5f6e78ca4701cac90131a0a718b00a3552",
    "type_VII_T_exp.fas": "e33d855e8d210fc3a29206441c9e7ce0eed147c830801941dcf6bfc531c0e7a1",
    "type_VIII_AT_exp.fas": "70373fe422199ce03b08f1866b95802944060a09f17c1b47213c1010862227d9",
    "type_VIII_T_exp.fas": "b85d3a26972e72b1ff44d4631b27f4d3b07c359d808b3d7e43a320874030ecfb",
}
ORITDB_URL = "https://bioinfo-mml.sjtu.edu.cn/oriTDB2/download/sequence/"
ORITDB_FILES = {
    "relaxase_all.fasta": "208e6f6d2e7376608d16b75d913f72bcc254f7d5cd2fe162c3946a038b23a91c",
    "auxiliary_all.fasta": "c6f2ceb0ab0b97e6348fc67b40e4d19648720d05545b4ba9c7d3f4fcff2ab974",
    "t4cp_all.fasta": "d36d05cf0f76d4b216f35cf44946268c0c1b8a5473f017b1026c9d791d1f45b0",
}
DBAPIS_URL = "https://pro.unl.edu/dbAPIS/download_file.php?file="  # bcb.unl.edu/dbAPIS redirects here
MOBILEOG_EVIDENCE = ("Manual", "Homology", "Keyword Search")  # header field 7
ALL_CLASSES = "all evidence classes kept, class recorded per label"

DBS = {
    "tadb": dict(
        version="TADB 3.0 (last update June 2023), experimentally validated protein entries",
        files=[(TADB_URL + n, n, h) for n, h in TADB_FILES.items()],
        url="https://bioinfo-mml.sjtu.edu.cn/TADB3/download.php",
        citation="Guan J, Chen Y, Goh YX, Wang M, Tai C, Deng Z, Song J, Ou HY. TADB 3.0: an "
                 "updated database of bacterial toxin-antitoxin loci and associated mobile "
                 "genetic elements. Nucleic Acids Res. 2024;52(D1):D784-D790. "
                 "doi:10.1093/nar/gkad962",
        licence="none stated on the website; the paper describes the database as freely "
                "available; article CC BY 4.0",
        evidence="experimentally validated entries only (the *_exp files)",
        header="<TADB id>|<file stem, e.g. type_II_T or regulator> <original TADB header>; "
               "the stem gives the TA type and the toxin (T) / antitoxin (AT) / regulator role"),
    "bacmet": dict(
        version="BacMet 2.0 (11 March 2018), BacMet2_EXP (753 experimentally confirmed genes)",
        files=[("http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta",
                "BacMet2_EXP_database.fasta",
                "e62756aafc1c653cdd1408fdbf558fd396ce1ef7ad99f2af50b0a6b67d00dd45"),
               ("http://bacmet.biomedicine.gu.se/download/BacMet2_EXP.753.mapping.txt",
                "BacMet2_EXP.753.mapping.txt",
                "d7f18a9a2b9bb5ae12a906dd37bef43147158e8ec43ecbe300f2c6776f374f8e")],
        url="http://bacmet.biomedicine.gu.se/download_temporary.html",
        citation="Pal C, Bengtsson-Palme J, Rensing C, Kristiansson E, Larsson DGJ. BacMet: "
                 "antibacterial biocide and metal resistance genes database. Nucleic Acids "
                 "Res. 2014;42(D1):D737-D743. doi:10.1093/nar/gkt1252 (no separate paper for "
                 "version 2.0)",
        licence="All rights reserved (website footer: Copyright 2013-2018 All rights "
                "reserved); no data licence stated; the paper describes the database as "
                "freely available",
        evidence="experimentally confirmed genes only (BacMet2_EXP)",
        header="original BacMet header, BacMet_ID|gene|...; compounds per BacMet_ID in "
               "raw/BacMet2_EXP.753.mapping.txt"),
    "oritdb": dict(
        version="oriTDB 2.0 (last update June 2024), relaxase, auxiliary protein and T4CP, "
                "experimentally validated plus predicted (_all)",
        files=[(ORITDB_URL + n, n, h) for n, h in ORITDB_FILES.items()],
        url="https://bioinfo-mml.sjtu.edu.cn/oriTDB2/download.php",
        citation="Liu G, Li X, Guan J, Tai C, Weng Y, Chen X, Ou HY. oriTDB: a database of the "
                 "origin-of-transfer regions of bacterial mobile genetic elements. Nucleic "
                 "Acids Res. 2025;53(D1):D163-D168. doi:10.1093/nar/gkae869",
        licence="none stated on the website; the paper describes the database as freely "
                "available; article CC BY 4.0",
        evidence="experimentally validated plus predicted entries (the _all files)",
        header="<role>_<n> <original oriTDB header>, role in relaxase, auxiliary, t4cp; "
               "the original header carries the MOB or T4CP family; oriTDB offers no T4SS "
               "protein download"),
    "card": dict(
        version=None,  # read from card.json
        files=[("https://card.mcmaster.ca/latest/data", "card-data.tar.bz2", None)],
        url="https://card.mcmaster.ca/latest/data",
        citation="Alcock BP, Huynh W, Chalil R, et al. CARD 2023: expanded curation, support "
                 "for machine learning, and resistome prediction at the Comprehensive "
                 "Antibiotic Resistance Database. Nucleic Acids Res. 2023;51(D1):D690-D699. "
                 "doi:10.1093/nar/gkac920",
        licence="free for non-commercial research or academic use by academic, government or "
                "non-profit institutions; commercial use requires a licence from McMaster "
                "University (card.mcmaster.ca/about, Terms of Use sections 4 and 5)",
        evidence="every protein homolog model",
        header="protein_fasta_protein_homolog_model.fasta as distributed by CARD; model "
               "cut-offs and ARO categories in card.json"),
    "mobileog": dict(
        version="mobileOG-db beatrix-1.6 (Zenodo 13241605, 2024-08-06), All set (Manual, "
                "Homology and Keyword Search entries)",
        files=[("https://zenodo.org/records/13241605/files/mobileOG-db_beatrix-1.6.All.faa",
                "mobileOG-db_beatrix-1.6.All.faa",
                "a6692f4d835355239643eb2592fe4c02c2a7824f9bcd14df3e0a50ec97dd74b4"),
               ("https://zenodo.org/records/13241605/files/mobileOG-db-beatrix-1.6.README.txt",
                "mobileOG-db-beatrix-1.6.README.txt",
                "145471e7a803878c51b1671708dce65a6627e90577e1a606e0cf626bc663f4c9")],
        url="https://zenodo.org/records/13241605 (archive deposited by the first author); "
            "website https://mobileogdb.flsi.cloud.vt.edu/ (Data Version: Beatrix 1.6 v1)",
        citation="Brown CL, Mullet J, Hindi F, Stoll JE, Gupta S, Choi M, Keenum I, "
                 "Vikesland P, Pruden A, Zhang L. mobileOG-db: a manually curated database "
                 "of protein families mediating the life cycle of bacterial mobile genetic "
                 "elements. Appl Environ Microbiol. 2022;88(18):e00991-22. "
                 "doi:10.1128/aem.00991-22",
        licence="Zenodo archive CC BY 4.0; GitHub repository GPL-3.0; article CC BY 4.0",
        evidence=ALL_CLASSES,
        header="original mobileOG header, id|name|UniProt|major category|minor "
               "categories|element|evidence (field 7: Manual, Homology or Keyword Search); "
               "headers can contain spaces (element 'Plasmid RefSeq', evidence 'Keyword "
               "Search'), so a search tool's subject id is the header up to the first space; "
               "119 entry ids occur twice, once as Homology and once as Keyword Search with "
               "the same sequence, and are told apart by that subject id"),
    "dbapis": dict(
        version="dbAPIS release 2026-06-18 (380 APIS families, 149 verified seeds), verified "
                "proteins and sequence homologues (anti_defense.pep)",
        files=[(DBAPIS_URL + "anti_defense.pep", "anti_defense.pep",
                "ee26941e609dcd25c8c862d4ed48474c0d9901d08be5a81e6377231643079320"),
               (DBAPIS_URL + "seed_and_familyrep_all_infor.tsv",
                "seed_and_familyrep_all_infor.tsv",
                "bd27f6096cd8c3a8e1b997eae84b10a65031924a3695f757304b18c0d435b6e6"),
               (DBAPIS_URL + "readme.txt", "readme.txt",
                "2fe51447889a611e01fdddbc9aa8e9c347d5e4ce7425c9976f84d2f1e90c7431")],
        url="https://bcb.unl.edu/dbAPIS (download page https://pro.unl.edu/dbAPIS/download.php; "
            "code at https://github.com/azureycy/dbAPIS, whose data_download holds the older "
            "2024-11-19 release)",
        citation="Yan Y, Zheng J, Zhang X, Yin Y. dbAPIS: a database of anti-prokaryotic "
                 "immune system genes. Nucleic Acids Res. 2024;52(D1):D419-D425. "
                 "doi:10.1093/nar/gkad932",
        licence="none stated (website footer: Copyright 2023 YIN LAB, UNL. All rights "
                "reserved; the GitHub repository has no licence file); article CC BY 4.0",
        evidence=ALL_CLASSES,
        header="<accession> gene=<gene> family=<family> evidence=verified|homolog. "
               "verified = the accession of a seed with if_verified = 1 in "
               "seed_and_familyrep_all_infor.tsv; gene = that seed's gene, or for a homologue "
               "the verified seed genes of its family (comma-separated; NA when the family "
               "has no verified seed); family = the APIS family, or the gene name of a "
               "singleton seed; 14 accessions assigned to two families with the same sequence "
               "are one record with both families comma-separated. The inhibited defence "
               "systems are in dbapis_metadata.tsv (protein_id, gene, family, "
               "defence_system, evidence), kept for reference and not read by the "
               "pipeline. The anti-CRISPR (Acr*) entries that "
               "anti_defense.pep also carries are not installed: dbAPIS states it does not "
               "cover anti-CRISPR proteins, and acrdb is the anti-CRISPR source"),
    "acrdb": dict(
        version="Anti-CRISPRdb v2.2 (March 2022), every entry of the core dataset "
                "(dataset1.csv)",
        files=[("https://web.archive.org/web/20240707082130id_/"
                "http://guolab.whu.edu.cn/anti-CRISPRdb/dataset/dataset1.csv",
                "dataset1.csv",
                "ba216c9f4401d60237755ac41496de6cdcfa7040dc6cdbf68951db1832a88c03")],
        url="http://cefg.uestc.cn/anti-CRISPRdb2 (host no longer resolves) and its successor "
            "http://guolab.whu.edu.cn/anti-CRISPRdb/ (HTTP 502 on 2026-09-25); the core "
            "dataset was taken from the Internet Archive capture of 2024-07-07; the file's "
            "SHA-1 equals the digest the archive recorded at capture",
        citation="Dong C, Wang X, Ma C, Zeng Z, Pu DK, Liu S, Wu CS, Chen S, Deng Z, Guo FB. "
                 "Anti-CRISPRdb v2.2: an online repository of anti-CRISPR proteins including "
                 "information on inhibitory mechanisms, activities and neighbors of curated "
                 "anti-CRISPR proteins. Database (Oxford). 2022;2022:baac010. "
                 "doi:10.1093/database/baac010",
        licence="none stated (website footer: Copyright CEFG 2021 All rights reserved); "
                "article CC BY 4.0",
        evidence=ALL_CLASSES,
        header="<anti_CRISPR_id> family=<Family> type=<Anti_type> acc=<Accession> "
               "evidence=<Verified|PLiterature|Putative>, the database's classify field"),
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(db, url, name, want, dest, local):
    """Place one raw file at dest, from a verified local copy or by download."""
    for d in local:
        cand = pathlib.Path(d) / db / name
        if want and cand.is_file() and sha256(cand) == want:
            shutil.copyfile(cand, dest)
            return f"local copy {cand} (SHA-256 verified)"
    part = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "plasmidann-installer"})
    with urllib.request.urlopen(req, timeout=600) as r, open(part, "wb") as out:
        shutil.copyfileobj(r, out, 1 << 20)
    got = sha256(part)
    if want and got != want:
        part.unlink()
        sys.exit(f"{db}/{name}: SHA-256 {got} does not match the pinned {want}; the file "
                 f"changed upstream at {url}; check it and update the pin")
    part.rename(dest)
    return f"download {url}"


def fasta(path):
    """(header, sequence) records; a line starting with whitespace before any sequence
    continues the header (oriTDB wraps some headers onto a second line)."""
    head, seq = None, []
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if head is not None:
                yield head, "".join(seq)
            head, seq = line[1:].strip(), []
        elif head is not None and not seq and line[:1].isspace():
            head += " " + line.strip()
        else:
            seq.append(line.strip())
    if head is not None:
        yield head, "".join(seq)


def write_fasta(path, records):
    """Write the records; an exact repeat of an earlier record (same header and sequence)
    is dropped, any other repeated identifier stops the install."""
    seen, n = {}, 0
    with open(path, "w") as out:
        for head, seq in records:
            seq = seq.rstrip("*").upper()
            sid = head.split()[0]
            digest = hashlib.sha256(f"{head}\n{seq}".encode()).digest()
            if sid in seen:
                if seen[sid] == digest:
                    continue
                sys.exit(f"{path}: identifier {sid} used for two different records")
            seen[sid] = digest
            out.write(f">{head}\n{seq}\n")
            n += 1
    if not n:
        sys.exit(f"{path}: no records")
    return n


def is_nucleotide(seq):
    return set(seq.upper()) <= set("ACGTUN")


def build_tadb(d, raw):
    dropped = 0

    def records():
        nonlocal dropped
        for name in TADB_FILES:
            stem = name[: -len("_exp.fas")]
            for head, seq in fasta(raw / name):
                if is_nucleotide(seq):  # RNA toxins and antitoxins of types I, III, VIII
                    dropped += 1
                    continue
                f = head.split(maxsplit=1)
                yield f"{f[0]}|{stem} {head}", seq

    n = write_fasta(d / "tadb.faa", records())
    return f"{n} protein entries; {dropped} nucleotide (RNA) entries excluded"


def build_bacmet(d, raw):
    n = write_fasta(d / "bacmet.faa", fasta(raw / "BacMet2_EXP_database.fasta"))
    return f"{n} entries"


def build_oritdb(d, raw):
    counts = {}

    def records():
        for name in ORITDB_FILES:
            role = name.split("_")[0]
            for i, (head, seq) in enumerate(fasta(raw / name), 1):
                counts[role] = i
                yield f"{role}_{i:05d} {head}", seq

    n = write_fasta(d / "oritdb.faa", records())
    return f"{n} entries ({', '.join(f'{k} {v}' for k, v in counts.items())})"


def build_card(d, raw):
    wanted = ("card.json", "protein_fasta_protein_homolog_model.fasta")
    with tarfile.open(raw / "card-data.tar.bz2") as tar:
        for m in tar.getmembers():
            if pathlib.PurePosixPath(m.name).name in wanted and m.isfile():
                with tar.extractfile(m) as src, open(d / pathlib.PurePosixPath(m.name).name,
                                                     "wb") as out:
                    shutil.copyfileobj(src, out)
    meta = json.load(open(d / "card.json"))
    n = sum(1 for _ in fasta(d / "protein_fasta_protein_homolog_model.fasta"))
    DBS["card"]["version"] = (f"CARD {meta['_version']} (card.json _timestamp "
                              f"{meta['_timestamp']}), protein homolog models")
    return f"{n} protein homolog model sequences"


def build_mobileog(d, raw):
    counts = dict.fromkeys(MOBILEOG_EVIDENCE, 0)

    def records():
        for head, seq in fasta(raw / "mobileOG-db_beatrix-1.6.All.faa"):
            f = head.split("|")
            if len(f) != 7 or f[6] not in counts:
                sys.exit(f"mobileOG header not id|name|uniprot|major|minor|element|"
                         f"evidence with a known evidence class: {head}")
            counts[f[6]] += 1
            yield head, seq

    n = write_fasta(d / "mobileog.faa", records())
    return f"{n} entries ({', '.join(f'{k} {v}' for k, v in counts.items())})"


def build_dbapis(d, raw):
    with open(raw / "seed_and_familyrep_all_infor.tsv", encoding="utf-8",
              errors="replace") as fh:
        table = list(csv.DictReader(fh, delimiter="\t"))
    verified = {r["Representative protein"]: r for r in table if r["if_verified"] == "1"}
    seed_genes, defence = {}, {}
    for r in table:
        fam = r["APIS families"]
        if r["if_verified"] == "1":
            seed_genes.setdefault(fam, []).append(r["APIS genes"])
        defence.setdefault(fam, [])
        if r["Defense systems"] and r["Defense systems"] not in defence[fam]:
            defence[fam].append(r["Defense systems"])

    # accession -> (families, sequence); an accession clustered into two families is one
    # record, and the anti-CRISPR entries (Acr*) are left to acrdb
    entries, acr = {}, 0
    for head, seq in fasta(raw / "anti_defense.pep"):
        fam, _, rest = head.split()[0].partition("|")
        if fam.startswith("Acr"):
            acr += 1
            continue
        seq = seq.rstrip("*").upper()
        fams, known = entries.setdefault(rest, ([], seq))
        if known != seq:
            sys.exit(f"dbAPIS accession {rest} carries two different sequences")
        if fam not in fams:
            fams.append(fam)
    missing = set(verified) - set(entries)
    if missing:
        sys.exit(f"dbAPIS verified seeds absent from anti_defense.pep: {sorted(missing)}")

    rows, counts = [], {"verified": 0, "homolog": 0}

    def records():
        for acc, (fams, seq) in entries.items():
            if acc in verified:
                evidence, genes = "verified", [verified[acc]["APIS genes"]]
            else:
                evidence = "homolog"
                genes = [g for f in fams for g in seed_genes.get(f, [])] or ["NA"]
            counts[evidence] += 1
            gene, family = ",".join(dict.fromkeys(genes)), ",".join(fams)
            rows.append((acc, gene, family,
                         "; ".join(x for f in fams for x in defence.get(f, [])), evidence))
            yield f"{acc} gene={gene} family={family} evidence={evidence}", seq

    n = write_fasta(d / "dbapis.faa", records())
    with open(d / "dbapis_metadata.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(("protein_id", "gene", "family", "defence_system", "evidence"))
        w.writerows(rows)
    return (f"{n} proteins (verified {counts['verified']}, homolog {counts['homolog']}) in "
            f"{len({f for r in rows for f in r[2].split(',')})} families; {acr} Acr entries "
            f"not installed")


def build_acrdb(d, raw):
    classes = {}

    def records():
        with open(raw / "dataset1.csv", encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh):
                classes[r["classify"]] = classes.get(r["classify"], 0) + 1
                # five sequences in dataset1.csv contain spaces, which DIAMOND rejects
                yield (f"{r['anti_CRISPR_id']} family={r['Family']} type={r['Anti_type']} "
                       f"acc={r['Accession']} evidence={r['classify']}",
                       "".join(r["Seq"].split()))

    n = write_fasta(d / "acrdb.faa", records())
    return f"{n} entries ({', '.join(f'{k} {v}' for k, v in classes.items())})"


BUILD = {"tadb": build_tadb, "bacmet": build_bacmet, "oritdb": build_oritdb,
         "card": build_card, "mobileog": build_mobileog, "dbapis": build_dbapis,
         "acrdb": build_acrdb}


def up_to_date(db, d):
    """True when VERSION, SOURCE and MANIFEST exist, the version is the expected one and
    every file in the manifest matches."""
    manifest = d / "MANIFEST.sha256"
    if not all((d / f).is_file() for f in ("VERSION", "SOURCE", "MANIFEST.sha256")):
        return False
    want = DBS[db]["version"]
    if want is not None and (d / "VERSION").read_text().strip() != want:
        return False
    for line in manifest.read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        p = d / name
        if not p.is_file() or sha256(p) != digest:
            return False
    return True


def install(db, labels, local):
    spec, d = DBS[db], labels / db
    if up_to_date(db, d):
        print(f"{db}: present and matches its manifest - skipped")
        return
    if d.exists():
        shutil.rmtree(d)
    raw = d / "raw"
    raw.mkdir(parents=True)
    origins = [fetch(db, url, name, want, raw / name, local)
               for url, name, want in spec["files"]]
    summary = BUILD[db](d, raw)
    (d / "VERSION").write_text(spec["version"] + "\n")
    (d / "SOURCE").write_text(
        f"url: {spec['url']}\n"
        f"citation: {spec['citation']}\n"
        f"licence: {spec['licence']}\n"
        f"evidence_rule: {spec['evidence']}\n"
        f"download_date: {datetime.date.today().isoformat()}\n"
        f"retrieved: {'; '.join(origins)}\n"
        f"content: {summary}\n"
        f"fasta_header: {spec['header']}\n"
        f"tool: tools/download_label_dbs.py\n")
    files = sorted(p for p in d.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256")
    (d / "MANIFEST.sha256").write_text(
        "".join(f"{sha256(p)}  {p.relative_to(d).as_posix()}\n" for p in files))
    print(f"{db}: installed - {summary}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", help="repository root")
    ap.add_argument("--local", action="append", default=[],
                    help="directory holding earlier downloads as DIR/<db>/<file name>; a "
                         "file is used only when its SHA-256 matches the pin (repeatable)")
    ap.add_argument("--db", action="append", choices=list(DBS),
                    help="install only these databases (default: all)")
    args = ap.parse_args()
    labels = pathlib.Path(args.root).resolve() / "data/refs/labels"
    for db in args.db or DBS:
        install(db, labels, args.local)


if __name__ == "__main__":
    main()
