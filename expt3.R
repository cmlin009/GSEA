# GSEA_pre-processing
library(openxlsx)
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
library(readr)

sig_genes_A20_vs_C1_deseq2_Result_5 <- read_csv("~/OneDrive/文件/丹蔘/RNAseq/Expressanalyst/A20_vs_C1/sig_genes_A20_vs_C1_deseq2_Result_5.csv")

GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'

#GSEA_KEGG
GSEA_input <- sig_genes_A20_vs_C1_deseq2_Result_5$logFC
names(GSEA_input) = sig_genes_A20_vs_C1_deseq2_Result_5$EntrezID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG, 'Are20_GSEA_KEGG.csv')

gseaplot2(GSEA_KEGG,8) #Adherens junction
gseaplot2(GSEA_KEGG,15) #Tight junction
gseaplot2(GSEA_KEGG,21) #Focal adhesion

#GSEA_GO
GSEA_GO <- gseGO(GSEA_input, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO, 'Are20_GSEA_GO.csv')

gseaplot2(GSEA_GO,6165) #epithelial to mesenchymal transition
gseaplot2(GSEA_GO,16) #cadherin binding
