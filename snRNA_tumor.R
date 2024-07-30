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

# CHRNA3 positive negative
View(Multiple.unpaired.t.tests.of.CHRNA3.tumor.mRNA.all.DEG)
CHRNA3 <- read.delim("~/OneDrive/文件/CHRNA/snRNAseq Puram 2017/GSEA/Multiple unpaired t tests of CHRNA3 tumor mRNA all DEG.txt")
CHRNA3 <- rename(CHRNA3, "X" = "Gene")
gene3 <- bitr(CHRNA3$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene3 <- rename (gene3, "SYMBOL" = 'Gene')
info3 <- left_join(CHRNA3, gene3, by = "Gene")

GSEA_input3 <- info3$Difference
names(GSEA_input3) = info3$ENTREZID
GSEA_input3 = sort(GSEA_input3, decreasing = TRUE)

GSEA_KEGG3 <- gseKEGG(GSEA_input3, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG3, 'CHRNA3_snRNA_GSEA_KEGG.csv')

GSEA_GO3 <- gseGO(GSEA_input3, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO3, 'CHRNA3_snRNA_GSEA_GO.csv')

gseaplot2(GSEA_KEGG3,244) #Natural killer cell mediated cytotoxicity
gseaplot2(GSEA_KEGG3,97) #PD-L1 expression and PD-1 checkpoint pathway in cancer
gseaplot2(GSEA_KEGG3,302) #Ras signaling pathway
gseaplot2(GSEA_KEGG3,296) #PI3K-Akt signaling pathway
gseaplot2(GSEA_GO3,1096) #focal adhesion assembly
gseaplot2(GSEA_GO3,4045) #cell junction assembly
gseaplot2(GSEA_GO3,6708) #epithelial to mesenchymal transition
gseaplot2(GSEA_GO3,3365) #adaptive immune response
gseaplot2(GSEA_GO3,7457) #positive regulation of cell cycle

# CHRNA5 positive negative
CHRNA5 <- read.delim("~/OneDrive/文件/CHRNA/snRNAseq Puram 2017/GSEA/Multiple unpaired t tests of CHRNA5 tumor mRNA all DEG.txt")
CHRNA5 <- rename(CHRNA5, "X" = "Gene")
gene5 <- bitr(CHRNA5$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene5 <- rename (gene5, "SYMBOL" = 'Gene')
info5 <- left_join(CHRNA5, gene5, by = "Gene")

GSEA_input5 <- info5$Difference
names(GSEA_input5) = info5$ENTREZID
GSEA_input5 = sort(GSEA_input5, decreasing = TRUE)

GSEA_KEGG5 <- gseKEGG(GSEA_input5, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG5, 'CHRNA5_snRNA_GSEA_KEGG.csv')

GSEA_GO5 <- gseGO(GSEA_input5, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO5, 'CHRNA5_snRNA_GSEA_GO.csv')

gseaplot2(GSEA_KEGG5,159) #Natural killer cell mediated cytotoxicity
gseaplot2(GSEA_KEGG5,179) #PD-L1 expression and PD-1 checkpoint pathway in cancer
gseaplot2(GSEA_KEGG5,116) #Ras signaling pathway
gseaplot2(GSEA_KEGG5,189) #PI3K-Akt signaling pathway
gseaplot2(GSEA_GO5,13) #focal adhesion
gseaplot2(GSEA_GO5,132) #cell-cell junction
gseaplot2(GSEA_GO5,2966) #epithelial to mesenchymal transition
gseaplot2(GSEA_GO5,4880) #adaptive immune response
gseaplot2(GSEA_GO5,334) #positive regulation of cell cycle

# CHRNA7 positive negative
CHRNA7 <- read.delim("~/OneDrive/文件/CHRNA/snRNAseq Puram 2017/GSEA/Multiple unpaired t tests of CHRNA7 tumor mRNA all DEG.txt")
CHRNA7 <- rename(CHRNA7, "X" = "Gene")
gene7 <- bitr(CHRNA7$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene7 <- rename (gene7, "SYMBOL" = 'Gene')
info7 <- left_join(CHRNA7, gene7, by = "Gene")

GSEA_input7 <- info7$Difference
names(GSEA_input7) = info7$ENTREZID
GSEA_input7 = sort(GSEA_input7, decreasing = TRUE)

GSEA_KEGG7 <- gseKEGG(GSEA_input7, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG7, 'CHRNA7_snRNA_GSEA_KEGG.csv')

GSEA_GO7 <- gseGO(GSEA_input7, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO7, 'CHRNA7_snRNA_GSEA_GO.csv')

gseaplot2(GSEA_KEGG7,120) #Natural killer cell mediated cytotoxicity
gseaplot2(GSEA_KEGG7,150) #PD-L1 expression and PD-1 checkpoint pathway in cancer
gseaplot2(GSEA_KEGG7,231) #Ras signaling pathway
gseaplot2(GSEA_KEGG7,98) #PI3K-Akt signaling pathway
gseaplot2(GSEA_GO7,3) #focal adhesion
gseaplot2(GSEA_GO7,76) #cell-cell junction
gseaplot2(GSEA_GO7,5049) #epithelial to mesenchymal transition
gseaplot2(GSEA_GO7,5375) #adaptive immune response
gseaplot2(GSEA_GO7,3004) #positive regulation of cell cycle

# CHRNA3 high low
CHRNA3 <- read.delim("~/OneDrive/文件/CHRNA/snRNAseq Puram 2017/GSEA_2_high_low/Multiple unpaired t tests of CHRNA3 tumor mRNA high low._刪除全0的基因.txt")
CHRNA3 <- rename(CHRNA3, "X" = "Gene")
gene3 <- bitr(CHRNA3$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene3 <- rename (gene3, "SYMBOL" = 'Gene')
info3 <- left_join(CHRNA3, gene3, by = "Gene")

GSEA_input3 <- info3$Difference
names(GSEA_input3) = info3$ENTREZID
GSEA_input3 = sort(GSEA_input3, decreasing = TRUE)

GSEA_KEGG3 <- gseKEGG(GSEA_input3, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG3, 'CHRNA3_snRNA_GSEA_KEGG_tumor_high_low.csv')

GSEA_GO3 <- gseGO(GSEA_input3, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO3, 'CHRNA3_snRNA_GSEA_GO_tumor_high_low.csv')

gseaplot2(GSEA_KEGG3,49) #Cell adhesion molecules
gseaplot2(GSEA_KEGG3,239) #Natural killer cell mediated cytotoxicity
gseaplot2(GSEA_KEGG3,272) #PD-L1 expression and PD-1 checkpoint pathway in cancer
gseaplot2(GSEA_KEGG3,11) #cell cycle
gseaplot2(GSEA_KEGG3,118) #Ras signaling pathway
gseaplot2(GSEA_KEGG3,113) #PI3K-Akt signaling pathway
gseaplot2(GSEA_GO3,5752) #cell-cell junction
gseaplot2(GSEA_GO3,6205) #epithelial to mesenchymal transition
gseaplot2(GSEA_GO3,2659) #adaptive immune response

# CHRNA5 high low
CHRNA5 <- read.delim("~/OneDrive/文件/CHRNA/snRNAseq Puram 2017/GSEA_2_high_low/Multiple unpaired t tests of CHRNA5 tumor mRNA high low._刪除全0的基因.txt")
CHRNA5 <- rename(CHRNA5, "X" = "Gene")
gene5 <- bitr(CHRNA5$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene5 <- rename (gene5, "SYMBOL" = 'Gene')
info5 <- left_join(CHRNA5, gene5, by = "Gene")

GSEA_input5 <- info5$Difference
names(GSEA_input5) = info5$ENTREZID
GSEA_input5 = sort(GSEA_input5, decreasing = TRUE)

GSEA_KEGG5 <- gseKEGG(GSEA_input5, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG5, 'CHRNA5_snRNA_GSEA_KEGG_tumor_high_low.csv')

GSEA_GO5 <- gseGO(GSEA_input5, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO5, 'CHRNA5_snRNA_GSEA_GO_tumor_high_low.csv')

gseaplot2(GSEA_KEGG5,72) #Cell adhesion molecules
gseaplot2(GSEA_KEGG5,83) #Natural killer cell mediated cytotoxicity
gseaplot2(GSEA_KEGG5,185) #PD-L1 expression and PD-1 checkpoint pathway in cancer
gseaplot2(GSEA_KEGG5,139) #cell cycle
gseaplot2(GSEA_KEGG5,183) #Ras signaling pathway
gseaplot2(GSEA_KEGG5,78) #PI3K-Akt signaling pathway
gseaplot2(GSEA_GO5,300) #cell-cell junction
gseaplot2(GSEA_GO5,7345) #epithelial to mesenchymal transition
gseaplot2(GSEA_GO5,233) #adaptive immune response

# CHRNA7 high low
CHRNA7 <- read.delim("~/OneDrive/文件/CHRNA/snRNAseq Puram 2017/GSEA_2_high_low/Multiple unpaired t tests of CHRNA7 tumor mRNA high low._刪除全0的基因.txt")
CHRNA7 <- rename(CHRNA7, "X" = "Gene")
gene7 <- bitr(CHRNA7$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene7 <- rename (gene7, "SYMBOL" = 'Gene')
info7 <- left_join(CHRNA7, gene7, by = "Gene")

GSEA_input7 <- info7$Difference
names(GSEA_input7) = info7$ENTREZID
GSEA_input7 = sort(GSEA_input7, decreasing = TRUE)

GSEA_KEGG7 <- gseKEGG(GSEA_input7, organism = KEGG_database, pvalueCutoff = 1)
write.csv(GSEA_KEGG7, 'CHRNA7_snRNA_GSEA_KEGG_tumor_high_low.csv')

GSEA_GO7 <- gseGO(GSEA_input7, ont = "ALL", OrgDb = GO_database, pvalueCutoff = 1)
write.csv(GSEA_GO7, 'CHRNA7_snRNA_GSEA_GO_tumor_high_low.csv')

gseaplot2(GSEA_KEGG7,118) #Cell adhesion molecules
gseaplot2(GSEA_KEGG7,69) #Natural killer cell mediated cytotoxicity
gseaplot2(GSEA_KEGG7,302) #PD-L1 expression and PD-1 checkpoint pathway in cancer
gseaplot2(GSEA_KEGG7,150) #cell cycle
gseaplot2(GSEA_GO7,1209) #cell-cell junction assembly
gseaplot2(GSEA_GO7,6466) #epithelial to mesenchymal transition
gseaplot2(GSEA_GO7,1797) #regulation of adaptive immune response
gseaplot2(GSEA_GO7,446) #protein kinase B signaling
gseaplot2(GSEA_GO7,7665) #regulation of small GTPase mediated signal transduction