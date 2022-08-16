library(seqinr)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)

setwd("/Volumes/lab-luscomben/working/slava/projects/tf-splicing")

humantfs = as_tibble(read.csv(file = "data/humantfs1.01/DatabaseExtract_v_1.01.csv", 
                              header = T, row.names = 1))
tfs = humantfs %>% 
  filter(Is.TF. == "Yes")

protein_fasta = read.fasta(file = "data/ensembl99/Homo_sapiens.GRCh38.pep.all.99.fa", seqtype = "AA")

is.valid.product = function(protein_record) {
  gene_id = str_match(getAnnot(protein_record), "gene:([A-Z0-9]+)")[2]
  transcript_biotype = str_match(getAnnot(protein_record), "transcript_biotype:([A-Za-z_]+)")[2]
  return((gene_id %in% tfs$Ensembl.ID) & (transcript_biotype == "protein_coding"))
}

protein_fasta_tfs = protein_fasta[unlist(lapply(protein_fasta, function(p) {is.valid.product(p)}))]

write.fasta(getSequence(protein_fasta_tfs), 
            str_sub(getAnnot(protein_fasta_tfs), start = 2), 
            "analysis/scan_interpro/output/protein_fasta/protein_fasta_tfs_99.fa")
