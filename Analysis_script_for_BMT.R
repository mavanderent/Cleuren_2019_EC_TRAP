### Ribo analysis
### Comparison of the Tie2_Cre with/without Bonemarrow transplant
### version 03_30_2019
rm(list=ls())

library(dplyr)
library(data.table)
library(stats)
library(preprocessCore)
library(gplots)
library(ggplot2)

### load expression data original_LPS
setwd("$PATH")
dat    <- read.csv("02_02_19_original_BMT_RPKM.csv", header = TRUE, row.names = 1)

### load metadata
setwd("$PATH")
info   <- read.csv("02_02_19_original_BMT_samples.csv", header = TRUE, row.names = 1)
info$mouse <- as.character(info$mouse)

### check validity
if (dim(dat)[2] != dim(info)[1]) {
  print ("Countmatrix and Samplesmatrix do not have same number of samples")
} else {
  if (sum(colnames(dat) == rownames(info)) == dim(dat)[2]) {
    print("Continue") 
  } else {
    print("Samplenames do not match")
  }
}


### 1st filter: Remove mitochondrial genes
toMatch1 <- "mt-"
noMT     <- dat[!(grepl(toMatch1,rownames(dat))),]
rm(dat,toMatch1)

### 2st filter: Remove Lars2 and Gm26924 genes because of problemetic expression lung BMT sample and SC samples
noMT     <- noMT[!(grepl("Lars2",rownames(noMT))),]
noMT     <- noMT[!(grepl("Gm26924",rownames(noMT))),]

### 3rd filter: Remove all with count over 1000000 (none present in this data set)
noMT     <- noMT[(apply(noMT, 1, max) < 1000000),]

### convert to TPM
exp     <- noMT[,rownames(info)]
samples <- info[colnames(noMT),]
n       <- dim(exp)[1]
p       <- dim(exp)[2]

for (i in 1:p) {exp[,i] = exp[,i] / sum(exp[,i]) * 1000000}

### filter out lowly expressed genes by group
info2     <- info[info$fraction=="bead", ]
exp2      <- exp[,rownames(info2) ] 

tissue    <- as.vector(info2$tissue)
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
rm(i,exp2,info2, IND, TPM, noMT, info, tissue, treatment, fraction, genotype, DF, DT, n, dat)
TPM_in           <- exp[,rownames(samples[samples$fraction=="input", ])] 
TPM_b            <- exp[,rownames(samples[samples$fraction=="bead", ])] 
colnames(TPM_in) <- paste("TPM",colnames(TPM_in),sep="_")
colnames(TPM_b)  <- paste("TPM",colnames(TPM_b),sep="_")

### Analysis
n1  <- dim(exp)[1]
exp <- log2(1+exp)

if (length(unique(samples$treatment))>1){
  L2_in           <- exp[,rownames(samples[samples$fraction=="input", ])] 
  L2_b            <- exp[,rownames(samples[samples$fraction=="bead", ])] 
  colnames(L2_in) <- paste("Log2",colnames(L2_in),sep="_")
  colnames(L2_b)  <- paste("Log2",colnames(L2_b),sep="_")
  samples_in       <- samples[samples$fraction=="input", ]
  samples_b        <- samples[samples$fraction=="bead", ]
}

### Analysis Ratio
de.samples = c()
de.genes = rep(F, n1)

for (genotype in unique(samples[,"genotype"])) {
  for (treatment in unique(samples[,"treatment"])) { 
    for (tissue in unique(samples[,"tissue"])) {
      sample.names = c()
      for (mouse in 1:max(samples[,"mouse"])) {
        sample.bead = rownames(samples)[samples[, "mouse"] == mouse & samples[, "genotype"] == genotype & samples[, "treatment"] == treatment & samples[, "tissue"] == tissue & samples[, "fraction"] == "bead"]
        sample.input = rownames(samples)[samples[, "mouse"] == mouse & samples[, "genotype"] == genotype & samples[, "treatment"] == treatment & samples[, "tissue"] == tissue & samples[, "fraction"] == "input"]
        if (length(sample.bead)==0 | length(sample.input)==0) next
        sample.name = paste(tissue,genotype,treatment, "de", mouse, sep = ".")
        exp[, sample.name] = exp[, sample.bead] -  exp[, sample.input]  
        sample.names = c(sample.names, sample.name)
        }
      de.samples = c(de.samples, sample.names)
      if (length(sample.names) > 1) {
      exp[, paste(tissue,genotype,treatment, "de", sep = ".")] = apply(exp[, sample.names, drop=F], 1, mean)  
      exp[, paste(tissue,genotype,treatment, "pvalue", sep = ".")] = 1.0
      for (i in 1:n1) {
        exp[i, paste(tissue,genotype,treatment, "pvalue", sep = ".")] = t.test(exp[i, sample.names])$p.value
        exp[is.na(exp[, paste(tissue,genotype,treatment, "pvalue", sep = ".")]), paste(tissue,genotype,treatment, "pvalue", sep = ".")] = 1.0
        exp[, paste(tissue,genotype,treatment, "qvalue", sep = ".")] = p.adjust(exp[, paste(tissue,genotype,treatment, "pvalue", sep = ".")], method="BH")
        exp[, paste(tissue,genotype,treatment, "de.genes", sep = ".")] = exp[, paste(tissue,genotype,treatment, "qvalue", sep = ".")] < 0.1 & abs(exp[, paste(tissue,genotype,treatment, "de", sep = ".")]) > 1
        de.genes = de.genes | exp[, paste(tissue,genotype,treatment, "de.genes", sep = ".")]
        }
      }
    }
  }
}
exp <-exp[,-(grep("de.genes",colnames(exp)))]
rm(i,mouse, sample.bead,sample.input,sample.name,tissue,treatment,sample.names,genotype)
ratio <- cbind(TPM_in,TPM_b, exp[,-(1:p)])
setwd("$PATH")
write.csv(ratio,"03_30_2019_analysis_BMT_ratio.csv")
rm(n1,p)

rm(de.samples,de.genes,exp)

if (length(unique(samples$treatment))>1){

### Analysis input
n1  <- dim(L2_in)[1]
p   <- dim(L2_in)[2]


de.samples <- c()
de.genes   <- rep(F, n1)
exp        <- L2_in
samples    <- samples_in

for (genotype in unique(samples[,"genotype"])) {
  for (tissue in unique(samples[,"tissue"])) {
    t.samples = c()
    t.grp = c()
    for (treatment in unique(samples[,"treatment"])) { 
      sample.names = c()
      grp.trt      = c()
      for (mouse in 1:max(samples[,"mouse"])) {
        sample.input = paste("Log2_",rownames(samples)[samples[, "mouse"] == mouse & samples[, "genotype"] == genotype & samples[, "treatment"] == treatment & samples[, "tissue"] == tissue & samples[, "fraction"] == "input"],sep="")
        grp.name     = samples$treatment[samples[, "mouse"] == mouse & samples[, "genotype"] == genotype & samples[, "treatment"] == treatment & samples[, "tissue"] == tissue & samples[, "fraction"] == "input"]
        if (length(sample.input)==0) next
        sample.name = paste(tissue,genotype,treatment, "de", mouse, sep = ".")
        exp[, sample.name] = exp[, sample.input]
        sample.names = c(sample.names, sample.name)
        grp.trt  = c(grp.trt,grp.name)
      }
      t.samples = c(t.samples, sample.names)
      t.grp     = c(t.grp, grp.trt)
      de.samples = c(de.samples, sample.names)
      if (length(sample.names) > 1) {
        exp[, paste(tissue,genotype,treatment, "de", sep = ".")] = apply(exp[, sample.names, drop=F], 1, mean) 
      }
    }
    exp[, paste(tissue,genotype,treatment, "L2FC", sep = ".")] = 1.0
    exp[, paste(tissue,genotype,treatment, "pvalue", sep = ".")] = 1.0
    
    for (i in 1:n1) {
    x <- data.frame(val = t(exp[i,t.samples]), grp = t.grp)
    y <- t.test(x[,1] ~ x$grp)
    exp[i, paste(tissue,genotype,treatment, "L2FC", sep = ".")] =as.numeric(y$estimate[2]-y$estimate[1])
    exp[i, paste(tissue,genotype,treatment, "pvalue", sep = ".")] = y$p.value
    exp[is.na(exp[, paste(tissue,genotype,treatment, "pvalue", sep = ".")]), paste(tissue,genotype,treatment, "pvalue", sep = ".")] = 1.0
    exp[, paste(tissue,genotype,treatment, "qvalue", sep = ".")] = p.adjust(exp[, paste(tissue,genotype,treatment, "pvalue", sep = ".")], method="BH")
    exp[, paste(tissue,genotype,treatment, "de.genes", sep = ".")] = exp[, paste(tissue,genotype,treatment, "qvalue", sep = ".")] < 0.1 & abs(exp[, paste(tissue,genotype,treatment, "L2FC", sep = ".")]) > 1
    de.genes = de.genes | exp[, paste(tissue,genotype,treatment, "de.genes", sep = ".")]
    }
  }
}

exp <-exp[,-(grep("de.genes",colnames(exp)))]
rm(i, x, y, mouse, sample.input,sample.name,tissue,treatment,sample.names,genotype,grp.name,grp.trt,t.grp,t.samples)
input <- cbind(TPM_in, exp[,-(1:p)])
setwd("$PATH")
write.csv(input,"03_30_2019_analysis_BMT_input.csv")
rm(n1,p)

rm(de.samples,de.genes,samples,exp,TPM_in,samples_in, L2_in)


### Analysis bead
n1  <- dim(L2_b)[1]
p   <- dim(L2_b)[2]


de.samples <- c()
de.genes   <- rep(F, n1)
exp        <- L2_b
samples    <- samples_b

for (genotype in unique(samples[,"genotype"])) {
  for (tissue in unique(samples[,"tissue"])) {
    t.samples = c()
    t.grp = c()
    for (treatment in unique(samples[,"treatment"])) { 
      sample.names = c()
      grp.trt      = c()
      for (mouse in 1:max(samples[,"mouse"])) {
        sample.bead = paste("Log2_",rownames(samples)[samples[, "mouse"] == mouse & samples[, "genotype"] == genotype & samples[, "treatment"] == treatment & samples[, "tissue"] == tissue & samples[, "fraction"] == "bead"],sep="")
        grp.name     = samples$treatment[samples[, "mouse"] == mouse & samples[, "genotype"] == genotype & samples[, "treatment"] == treatment & samples[, "tissue"] == tissue & samples[, "fraction"] == "bead"]
        if (length(sample.bead)==0) next
        sample.name = paste(tissue,genotype,treatment, "de", mouse, sep = ".")
        exp[, sample.name] = exp[, sample.bead]
        sample.names = c(sample.names, sample.name)
        grp.trt  = c(grp.trt,grp.name)
      }
      t.samples = c(t.samples, sample.names)
      t.grp     = c(t.grp, grp.trt)
      de.samples = c(de.samples, sample.names)
      if (length(sample.names) > 1) {
        exp[, paste(tissue,genotype,treatment, "de", sep = ".")] = apply(exp[, sample.names, drop=F], 1, mean) 
      }
    }
    exp[, paste(tissue,genotype,treatment, "L2FC", sep = ".")] = 1.0
    exp[, paste(tissue,genotype,treatment, "pvalue", sep = ".")] = 1.0
    for (i in 1:n1) {
      x <- data.frame(val = t(exp[i,t.samples]), grp = t.grp)
      y <- t.test(x[,1] ~ x$grp)
      exp[i, paste(tissue,genotype,treatment, "L2FC", sep = ".")] =as.numeric(y$estimate[2]-y$estimate[1])
      exp[i, paste(tissue,genotype,treatment, "pvalue", sep = ".")] = y$p.value
      exp[is.na(exp[, paste(tissue,genotype,treatment, "pvalue", sep = ".")]), paste(tissue,genotype,treatment, "pvalue", sep = ".")] = 1.0
      exp[, paste(tissue,genotype,treatment, "qvalue", sep = ".")] = p.adjust(exp[, paste(tissue,genotype,treatment, "pvalue", sep = ".")], method="BH")
      exp[, paste(tissue,genotype,treatment, "de.genes", sep = ".")] = exp[, paste(tissue,genotype,treatment, "qvalue", sep = ".")] < 0.1 & abs(exp[, paste(tissue,genotype,treatment, "L2FC", sep = ".")]) > 1
      de.genes = de.genes | exp[, paste(tissue,genotype,treatment, "de.genes", sep = ".")]
    }
  }
}

exp <-exp[,-(grep("de.genes",colnames(exp)))]
rm(i, x, y, mouse, sample.bead,sample.name,tissue,treatment,sample.names,genotype,grp.name,grp.trt,t.grp,t.samples)
bead <- cbind(TPM_b, exp[,-(1:p)])
setwd("$PATH")
write.csv(bead,"03_30_2019_analysis_BMT_bead.csv")
rm(n1,p)

rm(de.samples,de.genes,samples,exp,TPM_b,samples_b, L2_b)
}else{
  rm(samples,TPM_b,TPM_in)
}

