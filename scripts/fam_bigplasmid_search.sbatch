#!/bin/bash
# Does the most prevalent dark family (IMGPR_plasmid_3300022589_000018|2) also occur on the
# dark ORFs of BIG plasmids (>= 20 kb)? The dark-ORF extraction only ever covered the 71,414
# cryptic small plasmids, so the big-plasmid dark proteome has to be built first.
#
# Submit: sbatch scripts/fam_bigplasmid_search.sbatch
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 08:00:00
#SBATCH --mem=64G
#SBATCH -J fam_big
#SBATCH -o data/fam_3300022589/logs/fam_big_%j.out
set -euo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
export PATH=/gorilla/home/amujkic/.conda/envs/panaroo/bin:$PATH
D=data/fam_3300022589
PY=/gorilla/home/amujkic/.conda/envs/genesis/bin/python

# 1. Dark ORFs of every big plasmid, same extractor and same header grammar as the small set.
$PY scripts/extract_plasann_proteins.py --ids "$D/big_ids.txt" --out "$D/faa_big_dark" \
    --mode dark --jobs 24
cat "$D"/faa_big_dark/*.faa | awk '/^>/{keep = !seen[$0]++} keep' > "$D/big_dark_orfs.faa"
echo "big-plasmid dark ORFs: $(grep -c '^>' "$D/big_dark_orfs.faa")"

# 2. Search the family against them at the SAME criterion that defined the family in the first
#    place (30% identity, 80% bidirectional coverage), so "member" means one thing throughout.
mmseqs easy-search "$D/members.faa" "$D/big_dark_orfs.faa" "$D/hits_mmseqs.tsv" "$D/tmp_search" \
    --min-seq-id 0.3 -c 0.8 --cov-mode 0 -s 7.5 --max-seqs 5000 --threads 24 -v 1 \
    --format-output 'query,target,fident,alnlen,qcov,tcov,evalue,bits'
rm -rf "$D/tmp_search"
echo "mmseqs hits: $(wc -l < "$D/hits_mmseqs.tsv")"

# 3. Profile search as the sensitive counterpart: MMseqs at 30% id is a hard cut-off, an HMM
#    built from the family alignment will find diverged members the pairwise search misses.
hmmbuild --amino --cpu 24 -n fam3300022589 "$D/family.hmm" "$D/distinct_aln_i.faa"
hmmsearch --cpu 24 -E 1e-5 --domtblout "$D/hits_hmmer.domtbl" -o /dev/null \
    "$D/family.hmm" "$D/big_dark_orfs.faa"
echo "hmmer hits: $(grep -vc '^#' "$D/hits_hmmer.domtbl" || true)"
echo DONE
