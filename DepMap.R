# data import and pre-processing
library(tidyverse)
library(openxlsx)
library(readxl)
library(readr)

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

CHRNA3_depmap_sort_quatile <- read_excel("~/OneDrive/文件/CHRNA/depmap/GSEA/2nd/CHRNA3_depmap_sort_quatile.xlsx")
CHRNA3_depmap_sort_quatile <- as.data.frame(CHRNA3_depmap_sort_quatile)
CHRNA3_depmap_sort_quatile <- round_df(CHRNA3_depmap_sort_quatile, digits=0)
CHRNA3_depmap_sort_quatile <- rename (CHRNA3_depmap_sort_quatile, " " = "...1")
CHRNA3_depmap_sort_quatile_2 <- CHRNA3_depmap_sort_quatile[,-1]
rownames(CHRNA3_depmap_sort_quatile_2) <- CHRNA3_depmap_sort_quatile[,1]

CHRNA5_depmap_sort_quatile <- read_excel("~/OneDrive/文件/CHRNA/depmap/GSEA/2nd/CHRNA5_depmap_sort_quatile.xlsx")
CHRNA5_depmap_sort_quatile <- as.data.frame(CHRNA5_depmap_sort_quatile)
CHRNA5_depmap_sort_quatile <- round_df(CHRNA5_depmap_sort_quatile, digits=0)
CHRNA5_depmap_sort_quatile <- rename (CHRNA5_depmap_sort_quatile, " " = "...1")
CHRNA5_depmap_sort_quatile_2 <- CHRNA5_depmap_sort_quatile[,-1]
rownames(CHRNA5_depmap_sort_quatile_2) <- CHRNA5_depmap_sort_quatile[,1]

CHRNA7_depmap_sort_quatile <- read_excel("~/OneDrive/文件/CHRNA/depmap/GSEA/2nd/CHRNA7_depmap_sort_quatile.xlsx")
CHRNA7_depmap_sort_quatile <- as.data.frame(CHRNA7_depmap_sort_quatile)
CHRNA7_depmap_sort_quatile <- round_df(CHRNA7_depmap_sort_quatile, digits=0)
CHRNA7_depmap_sort_quatile <- rename (CHRNA7_depmap_sort_quatile, " " = "...1")
CHRNA7_depmap_sort_quatile_2 <- CHRNA7_depmap_sort_quatile[,-1]
rownames(CHRNA7_depmap_sort_quatile_2) <- CHRNA7_depmap_sort_quatile[,1]

# DESeq2
library(DESeq2)
coldata <- data.frame(condition = factor(rep(c('control', 'treat'), each = 15), levels = c('control', 'treat')))

dds <- DESeqDataSetFromMatrix(countData = CHRNA3_depmap_sort_quatile_2, colData = coldata, design= ~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds1 <- DESeq(dds)
res2 <- results(dds1, contrast = c('condition', 'treat', 'control'), independentFiltering=FALSE)
res2_2 <- data.frame(res2, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res2_2, 'deepmap_CHRNA3_quatile_DESeq2_nofilt.txt', col.names = NA, sep = '\t', quote = FALSE)

dds <- DESeqDataSetFromMatrix(countData = CHRNA5_depmap_sort_quatile_2, colData = coldata, design= ~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds1 <- DESeq(dds)
res2 <- results(dds1, contrast = c('condition', 'treat', 'control'), independentFiltering=FALSE)
res2_2 <- data.frame(res2, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res2_2, 'deepmap_CHRNA5_quatile_DESeq2_nofilt.txt', col.names = NA, sep = '\t', quote = FALSE)

dds <- DESeqDataSetFromMatrix(countData = CHRNA7_depmap_sort_quatile_2, colData = coldata, design= ~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds1 <- DESeq(dds)
res2 <- results(dds1, contrast = c('condition', 'treat', 'control'), independentFiltering=FALSE)
res2_2 <- data.frame(res2, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res2_2, 'deepmap_CHRNA7_quatile_DESeq2_nofilt.txt', col.names = NA, sep = '\t', quote = FALSE)

# GSEA_pre-processing
library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(GOplot)
library(DOSE)
library(ggnewscale)
library(topGO)
library(circlize)
library(ComplexHeatmap)
library(ggridges)

GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'

deepmap_CHRNA3_quatile_DESeq2_nofilt <- read_excel("~/OneDrive/文件/CHRNA/depmap/GSEA/2nd/CHRNA3的差異基因/deepmap_CHRNA3_quatile_DESeq2_nofilt.xlsx")
deepmap_CHRNA3_quatile_DESeq2_nofilt <- rename (deepmap_CHRNA3_quatile_DESeq2_nofilt, "Gene" = "...1")
gene <- bitr(deepmap_CHRNA3_quatile_DESeq2_nofilt$Gene,fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb = GO_database)
gene <- rename (gene, 'ENSEMBL' = "Gene")
deepmap_CHRNA3 <- left_join(deepmap_CHRNA3_quatile_DESeq2_nofilt, gene)

deepmap_CHRNA5_quatile_DESeq2_nofilt <- read_excel("~/OneDrive/文件/CHRNA/depmap/GSEA/2nd/CHRNA5的差異基因/deepmap_CHRNA5_quatile_DESeq2_nofilt.xlsx")
deepmap_CHRNA5_quatile_DESeq2_nofilt <- rename (deepmap_CHRNA5_quatile_DESeq2_nofilt, "Gene" = "...1")
gene <- bitr(deepmap_CHRNA5_quatile_DESeq2_nofilt$Gene,fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb = GO_database)
gene <- rename (gene, 'ENSEMBL' = "Gene")
deepmap_CHRNA5 <- left_join(deepmap_CHRNA5_quatile_DESeq2_nofilt, gene)

deepmap_CHRNA7_quatile_DESeq2_nofilt <- read_excel("~/OneDrive/文件/CHRNA/depmap/GSEA/2nd/CHRNA7的差異基因/deepmap_CHRNA7_quatile_DESeq2_nofilt.xlsx")
deepmap_CHRNA7_quatile_DESeq2_nofilt <- rename (deepmap_CHRNA7_quatile_DESeq2_nofilt, "...1" = "Gene")
gene <- bitr(deepmap_CHRNA7_quatile_DESeq2_nofilt$Gene,fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb = GO_database)
gene <- rename (gene, 'ENSEMBL' = "Gene")
deepmap_CHRNA3 <- left_join(deepmap_CHRNA7_quatile_DESeq2_nofilt, gene)

# GSEA_KEGG_GO_CHRNA3
GSEA_input <- deepmap_CHRNA3$log2FoldChange
names(GSEA_input) = deepmap_CHRNA3$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)

GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG, 'depmap_CHRNA3_KEGG_enrich_quatile.csv')

GSEA_GO <- gseGO(GSEA_input, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO, 'depmap_CHRNA3_GO_enrich_quatile.csv')

gseaplot2(GSEA_KEGG,24)
gseaplot2(GSEA_KEGG,85)
gseaplot2(GSEA_KEGG,135)
gseaplot2(GSEA_KEGG,247)
gseaplot2(GSEA_KEGG,133)
gseaplot2(GSEA_KEGG,216)
gseaplot2(GSEA_GO,1382)
gseaplot2(GSEA_GO,7807)
gseaplot2(GSEA_GO,7565)

# GSEA_KEGG_GO_CHRNA5
GSEA_input <- deepmap_CHRNA5$log2FoldChange
names(GSEA_input) = deepmap_CHRNA5$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)

GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG, 'depmap_CHRNA5_KEGG_enrich_quatile.csv')

GSEA_GO <- gseGO(GSEA_input, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO, 'depmap_CHRNA5_GO_enrich_quatile.csv')

gseaplot2(GSEA_KEGG,89)
gseaplot2(GSEA_KEGG,39)
gseaplot2(GSEA_KEGG,256)
gseaplot2(GSEA_KEGG,1)
gseaplot2(GSEA_KEGG,267)
gseaplot2(GSEA_KEGG,201)
gseaplot2(GSEA_GO,155)
gseaplot2(GSEA_GO,239)
gseaplot2(GSEA_GO,839)

# GSEA_KEGG_GO_CHRNA7
GSEA_input <- deepmap_CHRNA7$log2FoldChange
names(GSEA_input) = deepmap_CHRNA7$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)

GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG, 'depmap_CHRNA7_KEGG_enrich_quatile.csv')

GSEA_GO <- gseGO(GSEA_input, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO, 'depmap_CHRNA7_GO_enrich_quatile.csv')

gseaplot2(GSEA_KEGG,16)
gseaplot2(GSEA_KEGG,177)
gseaplot2(GSEA_KEGG,229)
gseaplot2(GSEA_KEGG,294)
gseaplot2(GSEA_KEGG,115)
gseaplot2(GSEA_KEGG,66)
gseaplot2(GSEA_GO,81)
gseaplot2(GSEA_GO,3393)
gseaplot2(GSEA_GO,5024)