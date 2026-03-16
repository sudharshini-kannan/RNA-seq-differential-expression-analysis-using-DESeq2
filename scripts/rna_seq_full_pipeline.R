# RNA-Seq Gene Expression Analysis Pipeline

# Dataset: GSE60450

# --------------------------------

# 1. Install and Load Libraries

# --------------------------------

if (!requireNamespace("BiocManager", quietly=TRUE))
install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly=TRUE))
BiocManager::install("DESeq2")

if (!requireNamespace("R.utils", quietly=TRUE))
install.packages("R.utils")

library(DESeq2)
library(R.utils)

# --------------------------------

# 2. Create Project Structure

# --------------------------------

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# --------------------------------

# 3. Download RNA-seq Dataset

# --------------------------------

download.file(
"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/suppl/GSE60450_Lactation-GenewiseCounts.txt.gz",
destfile = "data/GSE60450_counts.txt.gz"
)

# Unzip file

gunzip("data/GSE60450_counts.txt.gz", overwrite = TRUE)

# --------------------------------

# 4. Load Count Matrix

# --------------------------------

counts <- read.table(
"data/GSE60450_counts.txt",
header = TRUE,
sep = "\t",
row.names = 1
)

# Remove gene length column

counts <- counts[, -1]

# --------------------------------

# 5. Prepare Metadata

# --------------------------------

samples <- colnames(counts)

condition <- ifelse(grepl("-D", samples),
"pregnancy",
"lactation")

metadata <- data.frame(
row.names = samples,
condition = condition
)

# Fix column name formatting

colnames(counts) <- gsub("\.", "-", colnames(counts))

# --------------------------------

# 6. Create DESeq2 Dataset

# --------------------------------

dds <- DESeqDataSetFromMatrix(
countData = counts,
colData = metadata,
design = ~ condition
)

# --------------------------------

# 7. Differential Expression

# --------------------------------

dds <- DESeq(dds)

res <- results(dds)

write.csv(as.data.frame(res),
"results/deseq2_results.csv")

# --------------------------------

# 8. PCA Plot

# --------------------------------

vsd <- vst(dds, blind = FALSE)

pdf("results/PCA_plot.pdf")
plotPCA(vsd, intgroup = "condition")
dev.off()

# --------------------------------

# 9. Volcano Plot

# --------------------------------

res_df <- as.data.frame(res)

pdf("results/volcano_plot.pdf")

plot(res_df$log2FoldChange,
-log10(res_df$pvalue),
pch = 20,
xlab = "Log2 Fold Change",
ylab = "-Log10 p-value",
main = "Volcano Plot")

dev.off()

# --------------------------------

# 10. Heatmap of Top Genes

# --------------------------------

resOrdered <- res[order(res$padj),]

topgenes <- rownames(resOrdered)[1:20]

mat <- assay(vsd)[topgenes,]

pdf("results/heatmap_top_genes.pdf")

heatmap(mat, scale = "row")

dev.off()

print("RNA-seq analysis completed successfully!")

