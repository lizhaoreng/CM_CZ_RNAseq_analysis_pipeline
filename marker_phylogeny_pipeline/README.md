# Multilocus marker extraction and phylogeny pipeline

This pipeline was used to extract and analyze multilocus phylogenetic markers from in vitro RNA-seq Trinity assemblies and reference genomes of two *Cercospora* isolates.

## Markers

The following loci were used:

- ITS
- TEF1
- ACT
- CAL
- HIS3

## Main steps

1. Merge in vitro RNA-seq reads and assemble transcriptomes with Trinity.
2. Download reference marker sequences from NCBI using Entrez Direct.
3. Select curated reference marker sequences.
4. Extract marker sequences from Trinity assemblies.
5. Extract marker sequences from reference genomes.
6. Merge reference, Trinity-derived, and genome-derived markers.
7. Align each locus using MAFFT.
8. Trim alignments using trimAl.
9. Concatenate loci using AMAS.
10. Infer phylogeny using IQ-TREE.
11. Calculate pairwise p-distances.

## Requirements

- FastQC
- MultiQC
- fastp
- Trinity
- BLAST+
- SeqKit
- Entrez Direct
- MAFFT
- trimAl
- AMAS
- IQ-TREE
- Python 3

## Usage

Edit `config.sh` first.

```bash
bash run_all_marker_phylogeny.sh
```

Or run step by step:

```bash
bash 00_prepare_invitro_trinity.sh
bash 01_download_reference_markers.sh
bash 02_select_reference_sequences.sh
bash 03_extract_trinity_markers.sh
bash 04_extract_genome_markers.sh
bash 05_build_phylogeny_5gene.sh
python 06_calculate_pairwise_distance.py \
  concatenated/Cercospora_ITS_TEF1_ACT_CAL_HIS3_concat.fasta \
  pairwise_distance_5gene.tsv
```

## Main outputs

```text
iqtree_out/Cercospora_5gene_with_genomes.treefile
iqtree_out/Cercospora_5gene_with_genomes.contree
iqtree_noCAL/Cercospora_4locus_noCAL_with_genomes.treefile
iqtree_noCAL/Cercospora_4locus_noCAL_with_genomes.contree
pairwise_distance_5gene.tsv
```

## Notes

The extraction of Trinity-derived TEF1, ACT, CAL, and HIS3 marker regions was manually curated based on BLAST hits against reference marker sequences.
