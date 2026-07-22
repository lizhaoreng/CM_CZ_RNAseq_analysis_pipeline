#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Run the full marker extraction and multilocus phylogeny pipeline.
# ==============================================================================

echo "============================================================"
echo "Cercospora marker extraction and phylogeny pipeline"
echo "Start time: $(date)"
echo "============================================================"

bash 00_prepare_invitro_trinity.sh
bash 01_download_reference_markers.sh
bash 02_select_reference_sequences.sh
bash 03_extract_trinity_markers.sh
bash 04_extract_genome_markers.sh
bash 05_build_phylogeny_5gene.sh

python 06_calculate_pairwise_distance.py \
    concatenated/Cercospora_ITS_TEF1_ACT_CAL_HIS3_concat.fasta \
    pairwise_distance_5gene.tsv

echo "============================================================"
echo "Pipeline completed"
echo "End time: $(date)"
echo "============================================================"

echo "Main outputs:"
echo "  iqtree_out/Cercospora_5gene_with_genomes.treefile"
echo "  iqtree_out/Cercospora_5gene_with_genomes.contree"
echo "  iqtree_noCAL/Cercospora_4locus_noCAL_with_genomes.treefile"
echo "  iqtree_noCAL/Cercospora_4locus_noCAL_with_genomes.contree"
echo "  pairwise_distance_5gene.tsv"
