#!/usr/bin/env bash

#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 4G
#SBATCH --time 00:15:00

mkdir ../data/gtex8/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm_chunks

split --numeric-suffixes=1 \
    --suffix-length=3 \
    --lines=1000 \
    ../data/gtex8/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.tsv \
    ../data/gtex8/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm_chunks/chunk.
