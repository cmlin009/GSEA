# GSEA_pre-processing
library(tidyverse)
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

# GSEA_KEGG_GO_CHRNA3
dds <- read.table("./cbioportal_CHRNA3_DEG_quatile.tsv", sep = '\t', head=T)
gene <- bitr(dds$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene <- rename (gene, "SYMBOL" = 'Gene')
info <- left_join(dds, gene)

GSEA_input <- info$Log2.Ratio
names(GSEA_input) = info$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)

GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG, 'CHRNA3_KEGG_enrich_quatile.csv')

GSEA_GO <- gseGO(GSEA_input, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO, 'CHRNA3_GO_enrich_quatile.csv')

gseaplot2(GSEA_KEGG,14)
gseaplot2(GSEA_KEGG,281)
gseaplot2(GSEA_KEGG,232)
gseaplot2(GSEA_KEGG,314)
gseaplot2(GSEA_KEGG,75)
gseaplot2(GSEA_KEGG,92)
gseaplot2(GSEA_GO,1730)
gseaplot2(GSEA_GO,1826)
gseaplot2(GSEA_GO,2967)

# GSEA_KEGG_GO_CHRNA5
dds <- read.table("./cbioportal_CHRNA5_DEG_quatile.tsv", sep = '\t', head=T)
gene <- bitr(dds$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene <- rename (gene, "SYMBOL" = 'Gene')
info <- left_join(dds, gene)

GSEA_input <- info$Log2.Ratio
names(GSEA_input) = info$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)

GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG, 'CHRNA5_KEGG_enrich_quatile.csv')

GSEA_GO <- gseGO(GSEA_input, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO, 'CHRNA5_GO_enrich_quatile.csv')

gseaplot2(GSEA_KEGG,21)
gseaplot2(GSEA_KEGG,22)
gseaplot2(GSEA_KEGG,53)
gseaplot2(GSEA_KEGG,10)
gseaplot2(GSEA_KEGG,64)
gseaplot2(GSEA_KEGG,77)
gseaplot2(GSEA_GO,205)
gseaplot2(GSEA_GO,7184)
gseaplot2(GSEA_GO,13)

# GSEA_KEGG_GO_CHRNA7
dds <- read.table("./cbioportal_CHRNA7_DEG_quatile.tsv", sep = '\t', head=T)
gene <- bitr(dds$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene <- rename (gene, "SYMBOL" = 'Gene')
info <- left_join(dds, gene)

GSEA_input <- info$Log2.Ratio
names(GSEA_input) = info$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)

GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG, 'CHRNA7_KEGG_enrich_quatile.csv')

GSEA_GO <- gseGO(GSEA_input, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO, 'CHRNA7_GO_enrich_quatile.csv')

gseaplot2(GSEA_KEGG,3)
gseaplot2(GSEA_KEGG,183)
gseaplot2(GSEA_KEGG,150)
gseaplot2(GSEA_KEGG,197)
gseaplot2(GSEA_KEGG,58)
gseaplot2(GSEA_KEGG,17)
gseaplot2(GSEA_GO,1765)
gseaplot2(GSEA_GO,502)
gseaplot2(GSEA_GO,46)