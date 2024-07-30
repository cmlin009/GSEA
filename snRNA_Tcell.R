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

# CHRNA5 positive negative
CHRNA5_PN <- read.delim("~/OneDrive/文件/CHRNA/snRNAseq Puram 2017/GSEA_T cell/Multiple unpaired t tests of CHRNA5 T cell mRNA DEG positive negative 刪除全零.txt")
CHRNA5_PN <- rename(CHRNA5_PN, "X" = "Gene")
gene5_PN <- bitr(CHRNA5_PN$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene5_PN <- rename (gene5_PN, "SYMBOL" = 'Gene')
info5_PN <- left_join(CHRNA5_PN, gene5_PN, by = "Gene")

GSEA_input5_PN <- info5_PN$Difference
names(GSEA_input5_PN) = info5_PN$ENTREZID
GSEA_input5_PN = sort(GSEA_input5_PN, decreasing = TRUE)

GSEA_KEGG5_PN <- gseKEGG(GSEA_input5_PN, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG5_PN, 'CHRNA5_snRNA_GSEA_KEGG_Tcell_positive_negative.csv')

GSEA_GO5_PN <- gseGO(GSEA_input5_PN, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO5_PN, 'CHRNA5_snRNA_GSEA_GO_Tcell_positive_negative.csv')

gseaplot2(GSEA_KEGG5_PN,265) #PD-L1 expression and PD-1 checkpoint pathway in cancer
gseaplot2(GSEA_GO5_PN,3909) #adaptive immune response
gseaplot2(GSEA_GO5_PN,4602) #T cell activation

# CHRNA5 high low
CHRNA5_HL <- read.delim("~/OneDrive/文件/CHRNA/snRNAseq Puram 2017/GSEA_T cell/Multiple unpaired t tests of CHRNA5 T cell mRNA DEG high low 刪除全零.txt")
CHRNA5_HL <- rename(CHRNA5_HL, "X" = "Gene")
gene5_HL <- bitr(CHRNA5_HL$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene5_HL <- rename (gene5_HL, "SYMBOL" = 'Gene')
info5_HL <- left_join(CHRNA5_HL, gene5_HL, by = "Gene")

GSEA_input5_HL <- info5_HL$Difference
names(GSEA_input5_HL) = info5_HL$ENTREZID
GSEA_input5_HL = sort(GSEA_input5_HL, decreasing = TRUE)

GSEA_KEGG5_HL <- gseKEGG(GSEA_input5_HL, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG5_HL, 'CHRNA5_snRNA_GSEA_KEGG_Tcell_high_low.csv')

GSEA_GO5_HL <- gseGO(GSEA_input5_HL, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO5_HL, 'CHRNA5_snRNA_GSEA_GO_Tcell_high_low.csv')

gseaplot2(GSEA_KEGG5_HL,150) #PD-L1 expression and PD-1 checkpoint pathway in cancer
gseaplot2(GSEA_GO5_HL,380) #adaptive immune response
gseaplot2(GSEA_GO5_HL,77) #T cell activation

# CHRNA7 positive negative
CHRNA7_PN <- read.delim("~/OneDrive/文件/CHRNA/snRNAseq Puram 2017/GSEA_T cell/Multiple unpaired t tests of CHRNA7 T cell mRNA DEG positive negative 刪除全零.txt")
CHRNA7_PN <- rename(CHRNA7_PN, "X" = "Gene")
gene7_PN <- bitr(CHRNA7_PN$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene7_PN <- rename (gene7_PN, "SYMBOL" = 'Gene')
info7_PN <- left_join(CHRNA7_PN, gene7_PN, by = "Gene")

GSEA_input7_PN <- info7_PN$Difference
names(GSEA_input7_PN) = info7_PN$ENTREZID
GSEA_input7_PN = sort(GSEA_input7_PN, decreasing = TRUE)

GSEA_KEGG7_PN <- gseKEGG(GSEA_input7_PN, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG7_PN, 'CHRNA7_snRNA_GSEA_KEGG_Tcell_positive_negative.csv')

GSEA_GO7_PN <- gseGO(GSEA_input7_PN, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO7_PN, 'CHRNA7_snRNA_GSEA_GO_Tcell_positive_negative.csv')

gseaplot2(GSEA_KEGG7_PN,169) #PD-L1 expression and PD-1 checkpoint pathway in cancer
gseaplot2(GSEA_GO7_PN,1348) #adaptive immune response
gseaplot2(GSEA_GO7_PN,226) #T cell activation

# CHRNA7 high low
CHRNA7_HL <- read.delim("~/OneDrive/文件/CHRNA/snRNAseq Puram 2017/GSEA_T cell/Multiple unpaired t tests of CHRNA7 T cell mRNA DEG high low 刪除全零.txt")
CHRNA7_HL <- rename(CHRNA7_HL, "X" = "Gene")
gene7_HL <- bitr(CHRNA7_HL$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene7_HL <- rename (gene7_HL, "SYMBOL" = 'Gene')
info7_HL <- left_join(CHRNA7_HL, gene7_HL, by = "Gene")

GSEA_input7_HL <- info7_HL$Difference
names(GSEA_input7_HL) = info7_HL$ENTREZID
GSEA_input7_HL = sort(GSEA_input7_HL, decreasing = TRUE)

GSEA_KEGG7_HL <- gseKEGG(GSEA_input7_HL, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG7_HL, 'CHRNA7_snRNA_GSEA_KEGG_Tcell_high_low.csv')

GSEA_GO7_HL <- gseGO(GSEA_input7_HL, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO7_HL, 'CHRNA7_snRNA_GSEA_GO_Tcell_high_low.csv')

gseaplot2(GSEA_KEGG7_HL,60) #PD-L1 expression and PD-1 checkpoint pathway in cancer
gseaplot2(GSEA_GO7_HL,87) #adaptive immune response
gseaplot2(GSEA_GO7_HL,32) #T cell activation
