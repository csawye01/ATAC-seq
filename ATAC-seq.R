library("BiocInstaller")
source("http://bioconductor.org/biocLite.R")

# Load required packages
pkgNames <- c("rtracklayer", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene")
libSetup <- lapply(pkgNames, library, character.only = TRUE)

# Import narrowPeak file
fileName <- "atacseq_peaksnew.bed"
peakRanges <- import(fileName, format = "BED")

# Extract promoter ranges
promRanges <- promoters(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
as.data.frame(peakRanges)

# Map gene symbol names
promRanges$symbol <- mapIds(org.Hs.eg.db, promRanges$gene_id, "SYMBOL", "ENTREZID")

as.data.frame(promRanges)

# Find overlapping ranges
rangeOverlaps <- findOverlaps(peakRanges, promRanges)
peakRanges$symbol <- promRanges$symbol[subjectHits(rangeOverlaps)]
as.data.frame(rangeOverlaps)

ATACseq_peak_cnts <- read.table("atacseq_peak_counts.txt", sep="\t", row.names = "peakid",  header=TRUE)
ATACseq_design <- read.table("atacseq_design.txt", sep="\t", row.names = 1, header=TRUE)

# create DESeq
atacDDS <- DESeqDataSetFromMatrix(countData = ATACseq_peak_cnts, colData = ATACseq_design, design = ~condition)
atacDDS <- DESeq(atacDDS)
atac_Rlog <- rlog(atacDDS)

#Principle component analysis
plotPCA(atac_Rlog, intgroup = "condition", ntop = nrow(atac_Rlog))

res <- results(atacDDS)
res

#shortlist
inx1 <- res$padj < 0.05 & res$log2FoldChange > 1
inx2 <- res$padj < 0.05 & res$log2FoldChange < -1
res1 <- res[inx1,]
res2 <- res[inx2,]
res <- res[inx1 | inx2, ]
res <- res[complete.cases(res), ]

write.csv(res1, 
          file="ATAC-DESeq_results_padj0.05_lfc>1_upreg.csv")
write.csv(res2, 
          file="ATACDESeq_results_padj0.05_lfc<-1_downreg.csv")

upreg_highest <- read.table("ATAC-DESeq_1_upreg_chromloc.csv", sep=",",  header=FALSE)
downreg_highest <- read.table("ATACDESeq_1_downreg_chromreg.csv", sep=",", header=FALSE)

upreg.high.grange <- makeGRangesFromDataFrame(as.data.frame(upreg_highest),
                                              keep.extra.columns=FALSE,
                                              ignore.strand=FALSE,
                                              seqinfo=NULL,
                                              seqnames.field="V1",
                                              start.field="V2",
                                              end.field= "V3",
                                              strand.field="strand",
                                              starts.in.df.are.0based=FALSE)

downreg.high.grange <- makeGRangesFromDataFrame(as.data.frame(downreg_highest),
                                                keep.extra.columns=FALSE,
                                                ignore.strand=FALSE,
                                                seqinfo=NULL,
                                                seqnames.field="V1",
                                                start.field="V2",
                                                end.field= "V3",
                                                strand.field="strand",
                                                starts.in.df.are.0based=FALSE)
downreg.high.grange

downreg.overlaps <- subsetByOverlaps(promRanges, downreg.high.grange)
upreg.overlaps <- subsetByOverlaps(promRanges, upreg.high.grange)

library(biomaRt)
biocLite("VariantAnnotation")
library(VariantAnnotation)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno.upreg <- annotatePeak(upreg.overlaps, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno.upreg <- as.data.frame(peakAnno.upreg)
peakAnno.upreg$gene_id

peakAnno.downreg <- annotatePeak(downreg.overlaps, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno.downreg <- as.data.frame(peakAnno.downreg)


biocLite("clusterProfiler")
library(clusterProfiler)
## shows transcriptional programs from gene clustering

upreg.genes <- as.array(peakAnno.upreg$ENSEMBL) 
downreg.genes <- as.array(peakAnno.downreg$ENSEMBL)

ego <- enrichGO(comparison, 'org.Hs.eg.db', keyType = 'ENSEMBL', ont="BP")
dotplot(ego, showCategory=30)

ego.down <- enrichGO(downreg.genes, 'org.Hs.eg.db', keyType = 'ENSEMBL', ont="BP")
dotplot(ego.down, showCategory=30)


