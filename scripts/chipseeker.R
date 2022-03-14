BiocManager::install("ChIPseeker")
BiocManager::install(TxDb.Hsapiens.UCSC.hg38.knownGene)
BiocManager::install(EnsDb.Hsapiens.v86)
BiocManager::install(clusterProfiler)
BiocManager::install(AnnotationDbi)
BiocManager::install(org.Hs.eg.db)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

#loading consensus peaksets
#Change the file paths as per your requirements

files <- c("chip_seq_analysis/chipseq_r/consensus_peaks/H3K4me1_GM12878_overlaps_filt.broadPeak", "chip_seq_analysis/chipseq_r/consensus_peaks/Nfxl1_GM12878_idr0.05_filt.narrowPeak", "chip_seq_analysis/chipseq_r/consensus_peaks/Nfxl1_K562_idr0.05_filt.narrowPeak")
names(files) <- c("H3K4me1_GM12878", "Nfxl1_GM12878", "Nfxl1_K562")

#Getting Annotations for peaks
peakAnno <- lapply(files, annotatePeak, TxDb=edb,tssRegion=c(-2500, 2500), annoDb="org.Hs.eg.db", verbose=FALSE)
peakAnno

#barchart of genomic features for each sample
plotAnnoBar(peakAnno)

#Distribution of HistoneMark/TF-binding loci relative to TSS
plotDistToTSS(peakAnno, title="Distribution of histone mark or transcription factor-binding loci relative to TSS")

#Creating Dataframes of annotation
H3K4me1_GM12878_annot <- data.frame(peakAnno[["H3K4me1_GM12878"]]@anno)
head(H3K4me1_GM12878_annot)

Nfxl1_GM12878_annot <- data.frame(peakAnno[["Nfxl1_GM12878"]]@anno)
head(Nfxl1_GM12878_annot)

Nfxl1_K562_annot <- data.frame(peakAnno[["Nfxl1_K562"]]@anno)
head(Nfxl1_K562_annot)

#Getting EntreZ IDs for each sample
H3K4me1_GM12878_entrez <- H3K4me1_GM12878_annot$ENTREZID
Nfxl1_GM12878_entrez <- Nfxl1_GM12878_annot$ENTREZID
Nfxl1_K562_entrez <- Nfxl1_K562_annot$ENTREZID

#Running Gene-Ontology enrichment analysis
H3K4me1_GM12878_ego <- enrichGO(gene = H3K4me1_GM12878_entrez, 
                keyType = "ENTREZID", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05,
                readable = TRUE)

Nfxl1_GM12878_ego <- enrichGO(gene = Nfxl1_GM12878_entrez, 
                                keyType = "ENTREZID", 
                                OrgDb = org.Hs.eg.db, 
                                ont = "BP", 
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05, 
                                readable = TRUE)

Nfxl1_K562_ego <- enrichGO(gene = Nfxl1_K562_entrez, 
                                keyType = "ENTREZID", 
                                OrgDb = org.Hs.eg.db, 
                                ont = "BP", 
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05, 
                                readable = TRUE)

# Output results from GO analysis to a table
H3K4me1_GM12878_summary <- data.frame(H3K4me1_GM12878_ego)
write.csv(H3K4me1_GM12878_summary, "GOanalysis_H3K4me1_GM12878.csv")

Nfxl1_GM12878_summary <- data.frame(Nfxl1_GM12878_ego)
write.csv(Nfxl1_GM12878_summary, "GOanalysis_Nfxl1_GM12878.csv")

Nfxl1_K562_summary <- data.frame(Nfxl1_K562_ego)
write.csv(Nfxl1_K562_summary, "GOanalysis_Nfxl1_K562.csv")

# Dotplot visualization
dotplot(H3K4me1_GM12878_ego, showCategory=30, title="GO enrichment analysis of H3K4me1_GM12878")
dotplot(Nfxl1_GM12878_ego, showCategory=30, title="GO enrichment analysis of Nfxl1_GM12878")
dotplot(Nfxl1_K562_ego, showCategory=30, title="GO enrichment analysis of Nfxl1_K562")

#KEGG enrichment analysis
H3K4me1_GM12878_ekegg <- enrichKEGG(gene = H3K4me1_GM12878_entrez,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

Nfxl1_GM12878_ekegg <- enrichKEGG(gene = Nfxl1_GM12878_entrez,
                                    organism = 'hsa',
                                    pvalueCutoff = 0.05)

Nfxl1_K562_ekegg <- enrichKEGG(gene = Nfxl1_K562_entrez,
                                  organism = 'hsa',
                                  pvalueCutoff = 0.05)

# Dotplot visualization
dotplot(H3K4me1_GM12878_ekegg, showCategory=30, title="KEGG enrichment analysis of H3K4me1_GM12878")
dotplot(Nfxl1_GM12878_ekegg, showCategory=30, title="KEGG enrichment analysis of Nfxl1_GM12878")
dotplot(Nfxl1_K562_ekegg, showCategory=30, title="KEGG enrichment analysis of Nfxl1_K562")

#Comparing enrichment across samples
# Create a list with genes from each sample
genes_compare = lapply(c(peakAnno$Nfxl1_GM12878,peakAnno$Nfxl1_K562), function(i) as.data.frame(i)$ENTREZID)

# Run KEGG analysis
compKEGG <- compareCluster(geneCluster = list("Nfxl1_GM12878" = genes_compare[[1]], "Nfxl1_K562" = genes_compare[[2]]), 
                           fun = "enrichKEGG",
                           organism = "hsa",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis-Comparison")
