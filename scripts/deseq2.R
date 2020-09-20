# CardioPipeLine (c) 2020
# Author: Rohit Suratekar, ZDG Lab, IIMCB
#
# Functions related to DESeq2 analysis

suppressPackageStartupMessages(library("DESeq2"))

colData <- data.frame(sample = snakemake@params$samples,
                      condition = snakemake@params$conditions)

df <- read.csv(snakemake@input$combined)
rownames(df) <- df$gene_id
df$gene_id <- NULL

dds <- DESeq2::DESeqDataSetFromMatrix(countData = df,
                                      colData = colData,
                                      design = formula
                                      (snakemake@params$design))

# Relevel the condition by providing information about reference condition
print(dds$condition)
dds$condition <- relevel(dds$condition, ref = snakemake@params$reference)

print(dds$condition)

