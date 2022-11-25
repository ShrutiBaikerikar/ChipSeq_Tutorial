library(DiffBind)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)

register(SerialParam())

#load meta data of the samples
samples <- read.csv('D:/chipseq_r/meta_samples_diffbind.csv')
samples

#Reading the peaksets
nfxl.peaks <- dba(sampleSheet=samples)
nfxl.peaks
nfxl.peaks$config$doBlacklist = FALSE

#Applying greylist
bs_genome = BSgenome.Hsapiens.UCSC.hg38
GRCH38.ktype <- seqinfo(bs_genome)
nfxl.peaks <- dba.blacklist(nfxl.peaks, blacklist=FALSE, greylist = "BSgenome.Hsapiens.UCSC.hg38")


#Creating affinity binding matrix
nfxl.counts <- dba.count(nfxl.peaks)
nfxl.counts

#Normalizing
nfxl.counts <- dba.normalize(nfxl.counts)

#Differential Analysis
nfxl.model <- dba.contrast(nfxl.counts,design="~ Condition", contrast=c("Condition","Nfxl1_K562","Nfxl1_GM12878"), reorderMeta=list(Condition="Nfxl1_GM12878"))
nfxl.model

nfxl.model <- dba.analyze(nfxl.model)
dba.show(nfxl.model,bContrasts=TRUE)

nfxl.db <- dba.report(nfxl.model)
nfxl.db

#Plotting
dba.plotMA(nfxl.model)
dba.plotVolcano(nfxl.model)

