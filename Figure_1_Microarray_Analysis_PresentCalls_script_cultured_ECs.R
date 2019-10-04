### Micrarray PresentCalls file. 
### Last update 09_04_2019
### Changes: 
rm(list=ls())

library(lumi)
library(limma)
library(ggplot2)
library(dplyr)
library(reshape)

### LOAD DATA
setwd("$PATH")

fileName <- "03_03_2019_Raw_Data.csv"
dat   <- lumiR.batch(fileName)
rm(fileName)


### Add sample data 
meta <- read.csv("03_03_2019_sample_list.csv", header=TRUE, row.names = 1)
rownames(meta) == colnames(dat)
meta$group <- paste(meta$tissue,meta$day,sep="_")
meta$"sample" <- c(rep(c("A","B","C"),7),"A","A","A")
meta$"samplenr" <- rownames(meta)
rownames(meta) <- paste(meta$group,meta$sample,sep="_")
meta$sample <- NULL
colnames(dat) <- rownames(meta)
pData(dat) <- meta
rm(meta)

##  all samples, all genes, not normalized
dim(exprs(dat)) # 25697 genes    24 samples
### END OF LOAD DATA

## preprocessing
## all samples, all genes, expression conversion  # 25697 genes 24 samples

dat2 <- lumiExpresso(dat, bg.correct = TRUE, normalize = TRUE, variance.stabilize = TRUE, QC.evaluation = TRUE)
dim(exprs(dat2))
rm(dat)
### 

### PresentCall
dataKidney <- as.data.frame(exprs(dat2[,13:18]))
dataHeart  <- as.data.frame(exprs(dat2[,1:6]))

dat2_Kidney <- dat2[,13:18]
dat2_Heart  <- dat2[,1:6]
rm(dat2)

### To speed up the processing and reduce false positives, remove the unexpressed genes
presentCount_Kidney  <- as.data.frame(detectionCall(dat2_Kidney,type = "matrix"))
presentCount_Heart   <- as.data.frame(detectionCall(dat2_Heart, type = "matrix"))

## Add gene symbols to gene properties
probeList_Kidney <- fData(dat2_Kidney)
probeList_Heart  <- fData(dat2_Heart)

### Stich together
PresentCallList <- cbind(probeList_Kidney,presentCount_Kidney,presentCount_Heart,dataHeart,dataKidney)

rm(dat2_Heart,dat2_Kidney,dataHeart,dataKidney,presentCount_Heart,presentCount_Kidney,probeList_Heart,probeList_Kidney)
write.csv(PresentCallList, "PresentCalls_Microarray.csv")  
