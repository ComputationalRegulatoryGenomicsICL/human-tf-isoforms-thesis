#!/usr/bin/env bash

#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 32G
#SBATCH --time 12:00:00

fasta="../output/protein_fasta/protein_fasta_tfs_99.fa"

eval "$(conda shell.bash hook)"
conda activate ../../../tf-splicing/condaenv/tfsplicing/

../../../tools/interproscan/interproscan-5.39-77.0/interproscan.sh \
    --output-dir ../output/scan_results/all_analyses_ens99 \
    --disable-precalc \
    --input ${fasta} \
    --iprlookup \
    --tempdir ../temp

conda deactivate
