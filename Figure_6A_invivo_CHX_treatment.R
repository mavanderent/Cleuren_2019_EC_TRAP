### Principle component analysis; 
### with or without CHX (single cell protocol comparison) samples
### Version 09_03_2019

rm(list=ls())

### Load library
library(psych)
library(pcaMethods)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(SummarizedExperiment)
library(ggplot2)
library(cowplot)

### Load data
setwd("$PATH/")
### Log2 transformed TPMs
ratio       <- read.csv("09_03_2019_Ctrl_SC_bead_ECTx.csv", header = TRUE, row.names = 1)
### Sample info
meta_ratio  <- read.csv("09_03_2019_Ctrl_SC_samples_bead.csv", header = TRUE,row.names = 1)

sum(colnames(ratio)==rownames(meta_ratio))
meta_ratio$treatment <- as.factor( rep( c(rep("+ CHX",3),rep("- CHX",3)),3 )  )
meta_ratio$treatment <- factor(meta_ratio$treatment, levels = rev(levels(meta_ratio$treatment)))
meta_ratio$group     <- as.factor(paste(meta_ratio$tissue,meta_ratio$treatment,sep=" "))

### Plot
### MOCK RANGES
rowRanges <- GRanges(rep(c("chr1"), 14648),
                     IRanges(floor(runif(14648, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 14648, TRUE),
                     feature_id=sprintf("ID%03d", 1:14648))

ratio       <- SummarizedExperiment(assays=list(counts=as.matrix(ratio)), rowRanges = rowRanges, colData = DataFrame(meta_ratio))
ratio       <- DESeqTransform(ratio)
rm(rowRanges,meta_ratio)

rv <- rowVars(assay(ratio))
select <- order(rv, decreasing = TRUE)
pca <- prcomp(t(assay(ratio)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar_all <- round(percentVar[1:2]*100,digits = 0) 
my_x_lab <- paste0("PC1: ",percentVar_all[1], "% variance")
my_y_lab <- paste0("PC2: ",percentVar_all[2], "% variance")
rm(rv,select,pca,percentVar,percentVar_all)

data <-DESeq2::plotPCA(ratio, intgroup="group", returnData=TRUE)

pca_all <- qplot(PC1,PC2,shape=group ,color=group, data=data,size=I(2)) +
  scale_shape_manual(values=c(1,16,1,16,1,16)) +
  scale_color_manual(values = c("#E7298A","#E7298A", "#7570B3","#7570B3", "#1B9E77","#1B9E77")) +
  scale_y_continuous(limits = c(-50, 70)) + 
  scale_x_continuous(limits = c(-50, 70)) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(color="black", size=12 ),
        axis.title.y = element_text(color="black", size=12),
        legend.title =  element_blank()) +
  labs(x = my_x_lab, 
       y = my_y_lab) 

png("Figure_6A.png", width = 3, height = 3, units = "in", res = 600)
pca_all
dev.off()
rm(list=ls())