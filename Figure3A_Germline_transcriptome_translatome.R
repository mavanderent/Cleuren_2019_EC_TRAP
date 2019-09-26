### Principle component analysis 
### Germline input (transcriptome) vs bead (translatome)
### 09_03_2019
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

### Log2 transformed TPM
dat           <- read.csv("09_03_2019_germline_all.csv", header = TRUE, row.names = 1)
### samples
samples       <- read.csv("09_03_2019_samples_germline.csv", header = T,row.names = 1)
sum(colnames(dat)==rownames(samples))
samples$group <- paste(samples$tissue,samples$fraction,sep = " ")
samples       <- samples[,c(1,2,10,3:9)]

### Plot
### MOCK RANGES
rowRanges <- GRanges(rep(c("chr1"), 16631),
                     IRanges(floor(runif(16631, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 16631, TRUE),
                     feature_id=sprintf("ID%03d", 1:16631))

dat       <- SummarizedExperiment(assays=list(counts=as.matrix(dat)), rowRanges = rowRanges, colData = DataFrame(samples))
dat       <- DESeqTransform(dat)
rm(rowRanges)
rm(samples)

rv <- rowVars(assay(dat))
select <- order(rv, decreasing = TRUE)
pca <- prcomp(t(assay(dat)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar_all <- round(percentVar[1:2]*100,digits = 0) 
my_x_lab <- paste0("PC1: ",percentVar_all[1], "% variance")
my_y_lab <- paste0("PC2: ",percentVar_all[2], "% variance")
rm(rv,select,pca,percentVar,percentVar_all)

data <-DESeq2::plotPCA(dat, intgroup="group", returnData=TRUE)

pca <- qplot(PC1,PC2,shape=group, color=group, data=data, size=I(5)) +
  scale_shape_manual(values=rep(c(2,16),5)) +
  scale_color_manual(values = c("#E7298A","#E7298A", "#D95F02", "#D95F02", 
                                "#7570B3", "#7570B3","#1B9E77", "#1B9E77",
                                "#FFCC33","#FFCC33" )) +
  scale_y_continuous(limits = c((round(min(data$PC2))-10), (round(max(data$PC2))+10))) + 
  scale_x_continuous(limits = c((round(min(data$PC1))-10), (round(max(data$PC1))+10))) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(color="black", size=18 ),
        axis.title.y = element_text(color="black", size=18),
        legend.title =  element_blank()) +
  labs(x = my_x_lab, 
       y = my_y_lab)

png("Figure3A.png", width = 5, height = 4, units = "in", res = 600)
pca
dev.off()
rm(list=ls())
