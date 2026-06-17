library(dada2)
library(ggplot2)

path <- "/storage/praha1/home/solakh/data/papuang/novogene2/dem_trim_tail_18S3"  ## CHANGE ME to the directory containing the fastq files.
list.files(path)


# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))

#get sample names
sample.names <- basename(cutFs)
sample.names<-gsub("-trimmed_1.fastq.gz","",sample.names)
sample.names<-gsub("-assigned-","",sample.names)
tail(sample.names)

setwd("/storage/praha1/home/solakh/data/papuang/novogene2/dem_trim_tail_18S3")

#Names for filtered files
filtFs <- paste0(sample.names, "_READ1_filt.fastq.gz")
filtRs <- paste0(sample.names, "_READ2_filt.fastq.gz")
tail(filtRs)

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2,
    minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE) 
head(out)

write.csv(out, file = "filtering_summary.csv")

# Create OTU seqtab
fns <- list.files()
fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) 

fnFs <- fastqs[grepl("_READ1_filt.fastq.gz", fastqs)] 
fnRs <- fastqs[grepl("_READ2_filt.fastq.gz", fastqs)] 

sample.names <- gsub("_READ1_filt.fastq.gz","",fnFs)

#dereplicate sequences (get unique sequences from each file)
derepFs <- derepFastq(fnFs,n = 1e+05, verbose=T)
derepRs <- derepFastq(fnRs,n = 1e+05, verbose=T)

save(derepFs,file="/storage/praha1/home/solakh/data/papuang/novogene2/seqtab_18S_novo3/derepFs_PNG_18S_novo3.R")
save(derepRs,file="/storage/praha1/home/solakh/data/papuang/novogene2/seqtab_18S_novo3/derepRs_PNG_18S_novo3.R")

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#denoising
dadaFs <- dada(derepFs,err=NULL, multithread = T,selfConsist = TRUE,MAX_CONSIST=25)
dadaRs <- dada(derepRs,err=NULL, multithread = T,selfConsist = TRUE,MAX_CONSIST=25)

save(dadaFs,file="/storage/praha1/home/solakh/data/papuang/novogene2/seqtab_18S_novo3/dadaFs_18S_novo3.R")
save(dadaRs,file="/storage/praha1/home/solakh/data/papuang/novogene2/seqtab_18S_novo3/dadaRs_18S_novo3.R")

#merging
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = 10,maxMismatch=1,justConcatenate=F)

#final table
seqtab <- makeSequenceTable(mergers)

#save results
save(seqtab,file="/storage/praha1/home/solakh/data/papuang/novogene2/seqtab_18S_novo3/otutab_PNG_18S_novo3.R")
