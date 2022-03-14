BiocManager::install(“BiocParallel”)
BiocManager::install(“ChIPQC”)
BiocManager::install(“TxDb.Hsapiens.UCSC.hg38.knownGene”)

library(BiocParallel)
library(ChIPQC)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#load meta data of the samples
#Change the path to the file meta_samples_chipqc.csv to suit your system
samples <- read.csv('chip_seq_analysis/chipseq_r/meta_samples_chipqc.csv')
samples

#If you are using a Windows PC, please run this command
register(SerialParam())

## Create ChIPQC object [This step takes some time.]
chipObj <- ChIPQC(samples, annotation="hg38", consensus = TRUE, blacklist = "chip_seq_analysis/data/blacklist/GRCh38_unified_blacklist.bed" ) 

## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report", reportFolder="ChIPQCreport", facetBy = "Condition", lineBy = "Replicate")
