# GSEA_pre-processing
library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
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