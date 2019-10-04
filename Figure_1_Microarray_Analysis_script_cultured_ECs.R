### Micrarray analysis. 
### Last update 03_03_2019
### Changes: 
rm(list=ls())

library(lumi)
library(limma)
library(ggplot2)
library(dplyr)
library(reshape)

### LOAD DATA
setwd("$PATH/")

#Load expression file into lumi
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
# dim(exprs(dat)) # 25697 genes    24 samples
### END OF LOAD DATA

### START OF DIFFERENTIAL EXPRESSION ANALYSIS

## preprocessing
## all samples, all genes, expression conversion  # 21908 genes 24 samples

dat2 <- lumiExpresso(dat, bg.correct = TRUE, normalize = TRUE, variance.stabilize = TRUE, QC.evaluation = TRUE)
dat3 <- lumiExpresso(dat)

### 
rm(dat)
dataKidney <- exprs(dat2[,13:21])
dataHeart  <- exprs(dat2[,1:12])

dat2_Kidney <- dat2[,13:21]
dat2_Heart  <- dat2[,1:12]
rm(dat2)

### To speed up the processing and reduce false positives, remove the unexpressed genes
presentCount_Kidney  <- detectionCall(dat2_Kidney)
presentCount_Heart   <- detectionCall(dat2_Heart)

selDataKidney <- dataKidney[presentCount_Kidney > 0,]
selDataHeart  <- dataHeart[presentCount_Heart > 0,]

## Add gene symbols to gene properties
probeList_Kidney <- fData(dat2_Kidney)
probeList_Heart  <- fData(dat2_Heart)

probeList_Kidney <- probeList_Kidney[presentCount_Kidney > 0, ]
probeList_Heart  <- probeList_Heart[presentCount_Heart > 0, ]

rm(dat2_Heart,dat2_Kidney,dataHeart,dataKidney,presentCount_Heart,presentCount_Kidney)

## compare Kidney
sampleType <- c(rep("HD0",3), rep("HD3",3), rep("HD5",3))
design <- model.matrix(~ factor(sampleType) + 0 )
colnames(design) <- c("HD0","HD3", "HD5")
contrast.matrix<-makeContrasts(HD3.vs.HD0=HD3-HD0, HD5.vs.HD0=HD5-HD0, HD5.vs.HD3=HD5-HD3, levels=design)
fit <- lmFit(selDataKidney,design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

## Extract results
Kid_D3_D0 <- topTable(fit, coef="HD3.vs.HD0", adjust='BH', sort.by = "none", number=Inf)
Kid_D5_D0 <- topTable(fit, coef="HD5.vs.HD0", adjust='BH', sort.by = "none", number=Inf)
Kid_D5_D3 <- topTable(fit, coef="HD5.vs.HD3", adjust='BH', sort.by = "none", number=Inf)

colnames(Kid_D3_D0) <- paste("Kid_D3_D0",colnames(Kid_D3_D0),sep="_")
colnames(Kid_D5_D0) <- paste("Kid_D5_D0",colnames(Kid_D5_D0),sep="_")
colnames(Kid_D5_D3) <- paste("Kid_D5_D3",colnames(Kid_D5_D3),sep="_")

Kid_D3_D0 <- Kid_D3_D0[,c(1,2,4,5)]
Kid_D5_D0 <- Kid_D5_D0[,c(1,2,4,5)]
Kid_D5_D3 <- Kid_D5_D3[,c(1,2,4,5)]

Kidney <- cbind(probeList_Kidney,selDataKidney,Kid_D3_D0,Kid_D5_D0,Kid_D5_D3)
rm(Kid_D3_D0, Kid_D5_D0, Kid_D5_D3, probeList_Kidney,selDataKidney,design, contrast.matrix,fit, sampleType)
Kidney$ProbeID <- NULL
setwd("$PATH/")
write.csv(Kidney,"Differential_Expression_Kidney.csv")

## compare Heart
sampleType <- c(rep("HD0",3), rep("HD3",3), rep("HD5",3),rep("HD7",3))
design <- model.matrix(~ factor(sampleType)+0)
colnames(design) <- c("HD0","HD3", "HD5","HD7")
contrast.matrix<-makeContrasts(HD3.vs.HD0=HD3-HD0, HD5.vs.HD0=HD5-HD0, HD5.vs.HD3=HD5-HD3,
                               HD7.vs.HD0=HD7-HD0, HD7.vs.HD3=HD7-HD3, HD7.vs.HD5=HD7-HD5,
                               levels=design)
fit <- lmFit(selDataHeart,design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)


## Extract results
Heart_D3_D0 <- topTable(fit, coef="HD3.vs.HD0", adjust='BH', sort.by = "none", number=Inf)
Heart_D5_D0 <- topTable(fit, coef="HD5.vs.HD0", adjust='BH', sort.by = "none", number=Inf)
Heart_D5_D3 <- topTable(fit, coef="HD5.vs.HD3", adjust='BH', sort.by = "none", number=Inf)
Heart_D7_D0 <- topTable(fit, coef="HD7.vs.HD0", adjust='BH', sort.by = "none", number=Inf)
Heart_D7_D3 <- topTable(fit, coef="HD7.vs.HD3", adjust='BH', sort.by = "none", number=Inf)
Heart_D7_D5 <- topTable(fit, coef="HD7.vs.HD5", adjust='BH', sort.by = "none", number=Inf)

colnames(Heart_D3_D0) <- paste("Heart_D3_D0",colnames(Heart_D3_D0),sep="_")
colnames(Heart_D5_D0) <- paste("Heart_D5_D0",colnames(Heart_D5_D0),sep="_")
colnames(Heart_D5_D3) <- paste("Heart_D5_D3",colnames(Heart_D5_D3),sep="_")
colnames(Heart_D7_D0) <- paste("Heart_D7_D0",colnames(Heart_D7_D0),sep="_")
colnames(Heart_D7_D3) <- paste("Heart_D7_D3",colnames(Heart_D7_D3),sep="_")
colnames(Heart_D7_D5) <- paste("Heart_D7_D5",colnames(Heart_D7_D5),sep="_")

Heart_D3_D0 <- Heart_D3_D0[,c(1,2,4,5)]
Heart_D5_D0 <- Heart_D5_D0[,c(1,2,4,5)]
Heart_D5_D3 <- Heart_D5_D3[,c(1,2,4,5)]
Heart_D7_D0 <- Heart_D7_D0[,c(1,2,4,5)]
Heart_D7_D3 <- Heart_D7_D3[,c(1,2,4,5)]
Heart_D7_D5 <- Heart_D7_D5[,c(1,2,4,5)]

Heart <- cbind(probeList_Heart,selDataHeart,Heart_D3_D0,Heart_D5_D0,Heart_D5_D3, Heart_D7_D0, Heart_D7_D3, Heart_D7_D5)
rm(Heart_D7_D0, Heart_D7_D3, Heart_D7_D5, Heart_D3_D0, Heart_D5_D0, Heart_D5_D3, probeList_Heart,selDataHeart,design, contrast.matrix,fit, sampleType)
Heart$ProbeID <- NULL

write.csv(Heart,"Differential_Expression_Heart.csv")
rm(Heart, Kidney)
