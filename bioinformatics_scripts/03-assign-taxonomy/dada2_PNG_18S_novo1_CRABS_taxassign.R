library(ShortRead)
library(dada2)

rm(list = ls()) #clear the environment

#load the refseq
fasta<-readDNAStringSet("/storage/praha1/home/solakh/data/papuang/novogene2/seqtab_18S_novo1/REF_PNG_18S_novo1.fasta")

taxa = assignTaxonomy(as.character(fasta), "/storage/praha1/home/solakh/data/papuang/novogene2/seqtab_18S_novo1/CRABS/REF_dadaB_BlastFIlt.fasta" ,multithread=F,minBoot = 50)

save(taxa, file="/storage/praha1/home/solakh/data/papuang/novogene2/seqtab_18S_novo1/TAX_PNG_18S_novo1.R")
