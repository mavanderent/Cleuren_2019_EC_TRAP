### Ribo analysis HUVECS w/ Cycloheximide (CHX) & Harringtonine (HHT)
### version 09_03_2019
rm(list=ls())

library(dplyr)
library(data.table)
library(gplots)
library(ggplot2)
library(psych)
library(pcaMethods)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(SummarizedExperiment)
library(cowplot)

### load RPKM values data 
setwd("$PATH/")
dat    <- read.csv("09_02_19_RPKM_HUVEC_R1.csv", header = TRUE, row.names = 1)

### load metadata
info              <- read.csv("09_02_19_CHX_HHT_HUVEC_samples.csv", header = TRUE, row.names = 1)
info$mouse        <- as.character(rep(c(1,2,3,4),6))
colnames(info)[4] <- "human_sample"
info$treatment    <- as.character(c(rep("+ CHX",4),rep("+ HHT",4),rep(" Ctrl",4),rep("+ CHX",4),rep("+ HHT",4),rep(" Ctrl",4)))
info$cell         <- c(rep("HUVEC",12),rep("primary HUVEC (p2)",12))
info$group        <- as.factor( paste(info$cell,info$treatment,sep = " "))
info              <- info[,c(1:4,10,5,11,6:9)]

### check validity
if (dim(dat)[2] != dim(info)[1]) {
  print ("Countmatrix and Samplesmatrix do not same number of samples")
} else {
  if (sum(colnames(dat) == rownames(info)) == dim(dat)[2]) {
    print("Continue") 
  } else {
    print("Samplenames do not match")
  }
}


### 1st filter: Remove mitochondrial genes
toMatch1 <- "MT-"
noMT     <- dat[!(grepl(toMatch1,rownames(dat))),]
rm(dat,toMatch1)

### 2st filter: Remove Lars2 and Gm26924 genes because of problemetic expression lung BMT sample and SC samples
noMT     <- noMT[!(grepl("LARS2",rownames(noMT))),]

### 3rd filter: Remove all with count over 1000000 (none present in this data set)
noMT     <- noMT[(apply(noMT, 1, max) < 1000000),]

### convert to TPM
exp     <- noMT[,rownames(info)]
samples <- info[colnames(noMT),]
n       <- dim(exp)[1]
p       <- dim(exp)[2]

for (i in 1:p) {exp[,i] = exp[,i] / sum(exp[,i]) * 1000000}

### filter out lowly expressed genes by group
info2     <- info[info$fraction=="input", ]
exp2      <- exp[,rownames(info2) ] 

tissue    <- as.vector(info2$cell)
treatment <- as.vector(info2$treatment)
genotype  <- as.vector(info2$genotype)

TPM = 1

IND <- NULL

for (i in 1:dim(exp2)[1]){
  dat <- as.vector(t(exp2[i,]))
  DF  <- data.frame(tissue, genotype, treatment, dat)
  DT  <- data.table(DF)
  DT[,MEAN:=mean(dat), by=list(tissue, treatment, genotype)]
  if (sum(DT$MEAN>TPM)==0){
    IND[i] <- FALSE
  } else {
    IND[i] <- TRUE
  }
}

exp <- exp[IND,]
rm(i,exp2,info2, IND, TPM, noMT, info, tissue, treatment, genotype, DF, DT, n, dat,p)
TPM_in           <- exp[,rownames(samples[samples$fraction=="input", ])] 
colnames(TPM_in) <- paste("TPM",colnames(TPM_in),sep="_")
write.csv(TPM_in, "TPM_file.csv")
rm(TPM_in)

### Transformation
exp <- log2(1+exp)

### PCA plot

### separate different groups: HUVEC & primary HUVEC (p2)
sum(colnames(exp)==rownames(samples))
samples$treatment <- as.factor(samples$treatment)
samples$cell      <- as.factor(samples$cell)
samples$treatment <- ordered(samples$treatment, levels= c("Ctrl","CHX","HHT"))

exp_HUVEC         <- exp[,rownames(samples)]
samples_HUVEC     <- samples
rm(exp,samples)

### MOCK RANGES
rowRanges <- GRanges(rep(c("chr1"), 17141),
                     IRanges(floor(runif(17141, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 17141, TRUE),
                     feature_id=sprintf("ID%03d", 1:17141))

#### plot with input
HUVEC           <- SummarizedExperiment(assays=list(counts=as.matrix(exp_HUVEC)), rowRanges = rowRanges, colData = DataFrame(samples_HUVEC))
HUVEC           <- DESeqTransform(HUVEC)
rm(rowRanges,exp_HUVEC,samples_HUVEC)

rv_HUVECS <- rowVars(assay(HUVEC))
select_HUVECS <- order(rv_HUVECS, decreasing = TRUE)
pca_HUVEC <- prcomp(t(assay(HUVEC)[select_HUVECS, ]))
percentVar_HUVEC <- pca_HUVEC$sdev^2/sum(pca_HUVEC$sdev^2)
percentVar_HUVEC_all <- round(percentVar_HUVEC[1:2]*100,digits = 0) 

my_x_lab_HUVECS <- paste0("PC1: ",percentVar_HUVEC_all[1], "% variance")
my_y_lab_HUVECS <- paste0("PC2: ",percentVar_HUVEC_all[2], "% variance")

rm(rv_HUVECS,select_HUVECS,pca_HUVEC,percentVar_HUVEC,percentVar_HUVEC_all)

data_HUVEC <-DESeq2::plotPCA(HUVEC, intgroup="group", returnData=TRUE)

pca_HUVEC <- qplot(PC1, PC2, shape=group ,color=group, data=data_HUVEC,size=I(2)) +
  scale_shape_manual(values=c(16,1,15,16,1,15)) +
  scale_color_manual(values = c("#0033FF","#0033FF","#0033FF","#66CCFF","#66CCFF","#66CCFF")) +
  scale_y_continuous(limits = c(-40, 40)) + 
  scale_x_continuous(limits = c(-50, 50)) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(color="black", size=12 ),
        axis.title.y = element_text(color="black", size=12),
        legend.title =  element_blank()) +
  labs(x = my_x_lab_HUVECS, 
       y = my_y_lab_HUVECS) 

png("Figure_S1_.png", width = 4.5, height = 3, units = "in", res = 600)
pca_HUVEC
dev.off()
rm(list=ls())
