### Principle component analysis - Microarray; 
### Microarray 
rm(list=ls())
### Load library
library(psych)
library(pcaMethods)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(SummarizedExperiment)
library(ggplot2)

### Load data
setwd("$PATH/")

Kidney      <- read.csv("03_03_2019_Differential_Expression_Kidney.csv", header = TRUE, row.names = 1)
Heart       <- read.csv("03_03_2019_Differential_Expression_Heart.csv", header = TRUE, row.names = 1)
Kidney      <- Kidney[,2:10]
Heart       <- Heart[,2:13]
IND         <- rownames(Kidney)%in%rownames(Heart)
Kidney      <- Kidney[IND, ]
IND         <- rownames(Heart)%in%rownames(Kidney)
Heart       <- Heart[IND, ]
Kidney      <- Kidney[order(rownames(Kidney)), ]
Heart       <- Heart[order(rownames(Heart)), ]


input       <- cbind(Heart,Kidney)
rm(Kidney,Heart,IND)

meta        <- read.csv("03_03_2019_sample_list.csv", header = T,row.names = 1)
meta        <- meta[1:21,]
rownames(meta) <- colnames(input)

input_D0D3    <- input[,c(1:6,13:18)]
meta_D0D3     <- meta[c(1:6,13:18), ]
rm(input,meta)

### Plot
### MOCK RANGES
rowRanges   <- GRanges(rep("chr1", 12042),
                     IRanges(floor(runif(12042, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 12042, TRUE),
                     feature_id=sprintf("ID%03d", 1:12042))

input_D0D3  <- SummarizedExperiment(assays=list(counts=as.matrix(input_D0D3)), rowRanges = rowRanges, colData = DataFrame(meta_D0D3))
input_D0D3  <- DESeqTransform(input_D0D3)

rm(rowRanges)
rm(meta_D0D3)

png("Figure1.png", width = 5, height = 4, units = "in", res = 600)
pca_D0D3 <- plotPCA(input_D0D3, intgroup=c("tissue", "day"))
pca_D0D3 + geom_point(size=5,alpha=1) + 
  scale_color_manual(values = c("#D95F02","#F68B39","#7570B3","#A7A1FF"),
                     labels = c("Heart:Day0","Heart:Day3","Kidney:Day0","Kidney:Day3")) + 
  scale_y_continuous(limits = c(-30, 30)) + 
  scale_x_continuous(limits = c(-30,30)) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(color="black", size=18 ),
        axis.title.y = element_text(color="black", size=18))
dev.off()
rm(input_D0D3,pca_D0D3)

