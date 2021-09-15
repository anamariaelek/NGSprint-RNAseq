library(data.table)
library(stringr)
library(DESeq2)
library(msigdbr)
library(fgsea)

# counts Anastasia
fcounts_1 <- fread("featureCounts_results_Soroka.csv")
setnames(fcounts_1, c("gene", str_remove(colnames(fcounts_1)[-1],".bam")))
# counts Karol
fcounts_2 <- fread("FeatureCounts_all_reads.extracted.txt")
setnames(fcounts_2, c("gene", str_remove(colnames(fcounts_2)[-1],"(.subread)*.sorted.bam")))
# combine counts
fcounts <- merge.data.table(fcounts_1, fcounts_2, by="gene")

cts <- as.matrix(fcounts[,-1])
rownames(cts) <- fcounts$gene

# annotaion table
meta <- fread("../SRA_Run_Table.tsv")
meta_samples <- meta[match(colnames(fcounts)[-1],Run), .(sample,status,sex,age_sampling)]
coldata <- data.frame(
  meta_samples[,-1],
  row.names = meta_samples$sample
)
coldata$status <- factor(coldata$status, levels=c("CTRL","FRDA"))
coldata$sex <- factor(coldata$sex, levels=c("F","M"))

# DE analysis
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ sex + status
)
dds <- DESeq(dds)
res <- results(dds)

# shrink LFC values (for ranking and visualization)
resLFC <- lfcShrink(dds, coef=2, type="apeglm")
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

# viz MA plots
png("MA-plot.png", height=600, width=600, res=150)
par(mfrow=c(1,1))
plotMA(res, main="FRDA vs CTRL")
dev.off()

png("MA-plots-transform.png", height=800, width=800, res=150)
par(mfrow=c(2,2))
xlim <- c(1,1e5); ylim <- c(-4,4)
plotMA(res, xlim=xlim, ylim=ylim, main="no shrinkage")
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
dev.off()

# rank genes in the results
ranked_res <- res[order(res$padj),]
exampleRanks <- ranked_res$log2FoldChange
names(exampleRanks) <- rownames(ranked_res)
exampleRanks <- na.omit(exampleRanks) # reemove NAs!

# get gene sets
#tft_gene_sets = msigdbr(species = "Homo sapiens", category = "C3", subcategory = "GTRD")
cgp_gene_sets = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")
gene_sets <- unique(cgp_gene_sets$gs_name)
examplePathways <- lapply(gene_sets, function(x) cgp_gene_sets[cgp_gene_sets$gs_name==x,]$gene_symbol)
names(examplePathways) <- gene_sets

# GSEA
fgseaRes <- fgsea(
  pathways = examplePathways,
  stats = exampleRanks,
  minSize  = 15,
  maxSize  = 500
)

# plot GSEA results
pdf("FGSEA-C2-CP-pathways.pdf", height=8, width=10)
for (gs_name in gene_sets) {
  message(gs_name)
  examplePathway <- examplePathways[[gs_name]]
  if (length(examplePathways)>0)
    tryCatch(
      print(plotEnrichment(examplePathway, exampleRanks) + ggplot2::labs(title=gs_name)),
      error = function(e) warning(e)
    )
}
dev.off()

# table for top pathways
topPathwaysUp <- fgseaRes[ES > 0][padj<0.05][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][padj<0.05][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf("FGSEA-C2-CP-top-pathways-table.pdf", height=8, width=10)
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, gseaParam=1)
dev.off()
