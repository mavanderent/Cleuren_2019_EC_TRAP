### Version 09_03_2019
### Volcano plots
rm(list=ls())
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(grid)

setwd("$PATH/")

### load expression data
brain_input  <- read.csv("06_16_19_Ctrl_LPS_brain_input.csv", header = TRUE, row.names = 1)
brain_bead   <- read.csv("06_16_19_Ctrl_LPS_brain_ECTx.csv", header = TRUE, row.names = 1)
heart_input  <- read.csv("06_16_19_Ctrl_LPS_heart_input.csv", header = TRUE, row.names = 1)
heart_bead   <- read.csv("06_16_19_Ctrl_LPS_heart_ECTx.csv", header = TRUE, row.names = 1)
kidney_input <- read.csv("06_16_19_Ctrl_LPS_kidney_input.csv", header = TRUE, row.names = 1)
kidney_bead  <- read.csv("06_16_19_Ctrl_LPS_kidney_ECTx.csv", header = TRUE, row.names = 1)
liver_input  <- read.csv("06_16_19_Ctrl_LPS_liver_input.csv", header = TRUE, row.names = 1)
liver_bead   <- read.csv("06_16_19_Ctrl_LPS_liver_ECTx.csv", header = TRUE, row.names = 1)
lung_input   <- read.csv("06_16_19_Ctrl_LPS_lung_input.csv", header = TRUE, row.names = 1)
lung_bead    <- read.csv("06_16_19_Ctrl_LPS_lung_ECTx.csv", header = TRUE, row.names = 1)

databrain_Input <- data.frame(ensbl  = rownames(brain_input),
                              lfc    = brain_input$brain.input.L2FC,
                              padj   = -log10(brain_input$brain.input.qvalue))
databrain_Bead  <- data.frame(ensbl  = rownames(brain_bead),
                              lfc    = brain_bead$brain.bead.L2FC,
                              padj   = -log10(brain_bead$brain.bead.qvalue))
rm(brain_input)
rm(brain_bead)


dataheart_Input <- data.frame(ensbl  = rownames(heart_input),
                              lfc    = heart_input$heart.input.L2FC,
                              padj   = -log10(heart_input$heart.input.qvalue))
dataheart_Bead  <- data.frame(ensbl  = rownames(heart_bead),
                              lfc    = heart_bead$heart.bead.L2FC,
                              padj   = -log10(heart_bead$heart.bead.qvalue))
rm(heart_input)
rm(heart_bead)

datakidney_Input <- data.frame(ensbl  = rownames(kidney_input),
                              lfc    = kidney_input$kidney.input.L2FC,
                              padj   = -log10(kidney_input$kidney.input.qvalue))
datakidney_Bead  <- data.frame(ensbl  = rownames(kidney_bead),
                              lfc    = kidney_bead$kidney.bead.L2FC,
                              padj   = -log10(kidney_bead$kidney.bead.qvalue))
rm(kidney_input)
rm(kidney_bead)

dataliver_Input <- data.frame(ensbl  = rownames(liver_input),
                               lfc    = liver_input$liver.input.L2FC,
                               padj   = -log10(liver_input$liver.input.qvalue))
dataliver_Bead  <- data.frame(ensbl  = rownames(liver_bead),
                               lfc    = liver_bead$liver.bead.L2FC,
                               padj   = -log10(liver_bead$liver.bead.qvalue))
rm(liver_input)
rm(liver_bead)

datalung_Input <- data.frame(ensbl  = rownames(lung_input),
                              lfc    = lung_input$lung.input.L2FC,
                              padj   = -log10(lung_input$lung.input.qvalue))
datalung_Bead  <- data.frame(ensbl  = rownames(lung_bead),
                              lfc    = lung_bead$lung.bead.L2FC,
                              padj   = -log10(lung_bead$lung.bead.qvalue))
rm(lung_input)
rm(lung_bead)


# Modify dataset to add new colomn of colors
databrain_Input <- databrain_Input %>%
  mutate(color = ifelse((databrain_Input$padj>1 & databrain_Input$lfc > 1) | (databrain_Input$padj>1 & databrain_Input$lfc < -1.0), 
                        yes = "Differentially Expressed", 
                        no = "NS")) 

databrain_Bead <- databrain_Bead %>%
  mutate(color = ifelse((databrain_Bead$padj>1 & databrain_Bead$lfc > 1) | (databrain_Bead$padj>1 & databrain_Bead$lfc < -1.0), 
                        yes = "Differentially Expressed", 
                        no = "NS"))

dataheart_Input <- dataheart_Input %>%
  mutate(color = ifelse((dataheart_Input$padj>1 & dataheart_Input$lfc > 1) | (dataheart_Input$padj>1 & dataheart_Input$lfc < -1.0), 
                        yes = "Differentially Expressed", 
                        no = "NS")) 

dataheart_Bead <- dataheart_Bead %>%
  mutate(color = ifelse((dataheart_Bead$padj>1 & dataheart_Bead$lfc > 1) | (dataheart_Bead$padj>1 & dataheart_Bead$lfc < -1.0), 
                        yes = "Differentially Expressed", 
                        no = "NS"))

datakidney_Input <- datakidney_Input %>%
  mutate(color = ifelse((datakidney_Input$padj>1 & datakidney_Input$lfc > 1) | (datakidney_Input$padj>1 & datakidney_Input$lfc < -1.0), 
                        yes = "Differentially Expressed", 
                        no = "NS")) 

datakidney_Bead <- datakidney_Bead %>%
  mutate(color = ifelse((datakidney_Bead$padj>1 & datakidney_Bead$lfc > 1) | (datakidney_Bead$padj>1 & datakidney_Bead$lfc < -1.0), 
                        yes = "Differentially Expressed", 
                        no = "NS"))

dataliver_Input <- dataliver_Input %>%
  mutate(color = ifelse((dataliver_Input$padj>1 & dataliver_Input$lfc > 1) | (dataliver_Input$padj>1 & dataliver_Input$lfc < -1.0), 
                        yes = "Differentially Expressed", 
                        no = "NS")) 

dataliver_Bead <- dataliver_Bead %>%
  mutate(color = ifelse((dataliver_Bead$padj>1 & dataliver_Bead$lfc > 1) | (dataliver_Bead$padj>1 & dataliver_Bead$lfc < -1.0), 
                        yes = "Differentially Expressed", 
                        no = "NS"))

datalung_Input <- datalung_Input %>%
  mutate(color = ifelse((datalung_Input$padj>1 & datalung_Input$lfc > 1) | (datalung_Input$padj>1 & datalung_Input$lfc < -1.0), 
                        yes = "Differentially Expressed", 
                        no = "NS")) 

datalung_Bead <- datalung_Bead %>%
  mutate(color = ifelse((datalung_Bead$padj>1 & datalung_Bead$lfc > 1) | (datalung_Bead$padj>1 & datalung_Bead$lfc < -1.0), 
                        yes = "Differentially Expressed", 
                        no = "NS"))

max_lfc <- (round(max(c(databrain_Bead$lfc, databrain_Input$lfc, dataheart_Bead$lfc, dataheart_Input$lfc,
                        datakidney_Bead$lfc, datakidney_Input$lfc, dataliver_Bead$lfc, dataliver_Input$lfc,
                        datalung_Bead$lfc, datalung_Input$lfc))) + 1)

min_lfc <- (round(min(c(databrain_Bead$lfc, databrain_Input$lfc, dataheart_Bead$lfc, dataheart_Input$lfc,
                        datakidney_Bead$lfc, datakidney_Input$lfc, dataliver_Bead$lfc, dataliver_Input$lfc, 
                        datalung_Bead$lfc, datalung_Input$lfc))) - 1)

max_padj <- (round(max(c(databrain_Input$padj, databrain_Bead$padj, dataheart_Input$padj, dataheart_Bead$padj,
                         datakidney_Input$padj, datakidney_Bead$padj, dataliver_Input$padj, dataliver_Bead$padj,
                         datalung_Input$padj, datalung_Bead$padj))) + 1)


coloredBrain_Input <- ggplot(databrain_Input, aes(x = lfc, y = padj)) + 
  geom_point(aes(color = factor(color)), size = 1, alpha = 0.8, na.rm = T) + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5), legend.justification = c(0,1),
        legend.position = "none", panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour="black"),plot.margin = unit(c(0.1,0.1,0.2,0.1),"in")) + 
  xlab(expression(log[2]("LPS" / "control"))) + 
  ylab(expression(-log[10]("q-value"))) + 
  geom_hline(yintercept = 1., colour = "black", linetype="dashed") + 
  geom_vline(xintercept = -1.0, colour = "black", linetype="dashed") + 
  geom_vline(xintercept = 1.0, colour = "black", linetype="dashed") +
  scale_y_continuous(limits = c(0, max_padj)) + 
  scale_x_continuous(limits = c(min_lfc, max_lfc)) + 
  scale_color_manual(values = c("Differentially Expressed" = "#E64B35", 
                                "NS" = "#636363"),
                     labels = c("Differentially Expressed genes",
                                "Not significant")) 

coloredBrain_Bead <- ggplot(databrain_Bead, aes(x = lfc, y = padj)) + 
  geom_point(aes(color = factor(color)), size = 1, alpha = 0.8, na.rm = T) + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5), legend.justification = c(0,1),
        legend.position = "none", panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour="black"),plot.margin = unit(c(0.1,0.1,0.2,0.1),"in")) + 
  xlab(expression(log[2]("LPS" / "control"))) + 
  ylab(expression(-log[10]("q-value"))) + 
  
  geom_hline(yintercept = 1., colour = "black", linetype="dashed") + 
  geom_vline(xintercept = -1.0, colour = "black", linetype="dashed") + 
  geom_vline(xintercept = 1.0, colour = "black", linetype="dashed") +
  scale_y_continuous(limits = c(0, max_padj)) + 
  scale_x_continuous(limits = c(min_lfc, max_lfc)) + 
  scale_color_manual(values = c("Differentially Expressed" = "#E64B35", 
                                "NS" = "#636363"),
                     labels = c("Differentially Expressed genes",
                                "Not significant"))

coloredKidney_Input <- ggplot(datakidney_Input, aes(x = lfc, y = padj)) + 
  geom_point(aes(color = factor(color)), size = 1, alpha = 0.8, na.rm = T) + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5), legend.justification = c(0,1),
        legend.position = "none", panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour="black"),plot.margin = unit(c(0.1,0.1,0.2,0.1),"in")) + 
  xlab(expression(log[2]("LPS" / "control"))) +
  ylab(expression(-log[10]("q-value"))) +  
  geom_hline(yintercept = 1., colour = "black", linetype="dashed") + 
  geom_vline(xintercept = -1.0, colour = "black", linetype="dashed") + 
  geom_vline(xintercept = 1.0, colour = "black", linetype="dashed") +
  scale_y_continuous(limits = c(0, max_padj)) + 
  scale_x_continuous(limits = c(min_lfc, max_lfc)) + 
  scale_color_manual(values = c("Differentially Expressed" = "#E64B35", 
                                "NS" = "#636363"),
                     labels = c("Differentially Expressed genes",
                                "Not significant")) 

coloredKidney_Bead <- ggplot(datakidney_Bead, aes(x = lfc, y = padj)) + 
  geom_point(aes(color = factor(color)), size = 1, alpha = 0.8, na.rm = T) + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5), legend.justification = c(0,1),
        legend.position = "none", panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour="black"),plot.margin = unit(c(0.1,0.1,0.2,0.11),"in")) + 
  xlab(expression(log[2]("LPS" / "control"))) + 
  ylab(expression(-log[10]("q-value"))) + 
  geom_hline(yintercept = 1., colour = "black", linetype="dashed") + 
  geom_vline(xintercept = -1.0, colour = "black", linetype="dashed") + 
  geom_vline(xintercept = 1.0, colour = "black", linetype="dashed") +
  scale_y_continuous(limits = c(0, max_padj)) + 
  scale_x_continuous(limits = c(min_lfc, max_lfc)) + 
  scale_color_manual(values = c("Differentially Expressed" = "#E64B35", 
                                "NS" = "#636363"),
                     labels = c("Differentially Expressed genes",
                                "Not significant")) 

coloredLiver_Input <- ggplot(dataliver_Input, aes(x = lfc, y = padj)) + 
  geom_point(aes(color = factor(color)), size = 1, alpha = 0.8, na.rm = T) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5), legend.justification = c(0,1),
        legend.position = "none", panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour="black"),plot.margin = unit(c(0.1,0.1,0.2,0.1),"in")) + 
  xlab(expression(log[2]("LPS" / "control"))) + 
  ylab(expression(-log[10]("q-value"))) +  
  geom_hline(yintercept = 1., colour = "black", linetype="dashed") + 
  geom_vline(xintercept = -1.0, colour = "black", linetype="dashed") + 
  geom_vline(xintercept = 1.0, colour = "black", linetype="dashed") +
  scale_y_continuous(limits = c(0, max_padj)) + 
  scale_x_continuous(limits = c(min_lfc, max_lfc)) + 
  scale_color_manual(values = c("Differentially Expressed" = "#E64B35", 
                                "NS" = "#636363"),
                     labels = c("Differentially Expressed genes",
                                "Not significant")) 

coloredLiver_Bead <- ggplot(dataliver_Bead, aes(x = lfc, y = padj)) + 
  geom_point(aes(color = factor(color)), size = 1, alpha = 0.8, na.rm = T) + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5), legend.justification = c(0,1),
        legend.position = "none", panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour="black"),plot.margin = unit(c(0.1,0.1,0.2,0.1),"in")) + 
  xlab(expression(log[2]("LPS" / "control"))) + 
  ylab(expression(-log[10]("q-value"))) + 
  geom_hline(yintercept = 1., colour = "black", linetype="dashed") + 
  geom_vline(xintercept = -1.0, colour = "black", linetype="dashed") + 
  geom_vline(xintercept = 1.0, colour = "black", linetype="dashed") +
  scale_y_continuous(limits = c(0, max_padj)) + 
  scale_x_continuous(limits = c(min_lfc, max_lfc)) + 
  scale_color_manual(values = c("Differentially Expressed" = "#E64B35", 
                                "NS" = "#636363"),
                     labels = c("Differentially Expressed genes",
                                "Not significant")) 

coloredHeart_Input <- ggplot(dataheart_Input, aes(x = lfc, y = padj)) + 
  geom_point(aes(color = factor(color)), size = 1, alpha = 0.8, na.rm = T) + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5), legend.justification = c(0,1),
        legend.position = "none", panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour="black"),plot.margin = unit(c(0.1,0.1,0.2,0.1),"in")) + # clean up theme
  xlab(expression(log[2]("LPS" / "control"))) + 
  ylab(expression(-log[10]("q-value"))) + 
  geom_hline(yintercept = 1., colour = "black", linetype="dashed") + 
  geom_vline(xintercept = -1.0, colour = "black", linetype="dashed") + 
  geom_vline(xintercept = 1.0, colour = "black", linetype="dashed") +
  scale_y_continuous(limits = c(0, max_padj)) + 
  scale_x_continuous(limits = c(min_lfc, max_lfc)) +
  scale_color_manual(values = c("Differentially Expressed" = "#E64B35", 
                                "NS" = "#636363"),
                     labels = c("Differentially Expressed genes",
                                "Not significant")) 

coloredHeart_Bead <- ggplot(dataheart_Bead, aes(x = lfc, y = padj)) + 
  geom_point(aes(color = factor(color)), size = 1, alpha = 0.8, na.rm = T) + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5), legend.justification = c(0,1),
        legend.position = "none", panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour="black"),plot.margin = unit(c(0.1,0.1,0.2,0.1),"in")) + 
  xlab(expression(log[2]("LPS" / "control"))) + 
  ylab(expression(-log[10]("q-value"))) + 
  geom_hline(yintercept = 1., colour = "black", linetype="dashed") + 
  geom_vline(xintercept = -1.0, colour = "black", linetype="dashed") + 
  geom_vline(xintercept = 1.0, colour = "black", linetype="dashed") +
  scale_y_continuous(limits = c(0, max_padj)) + 
  scale_x_continuous(limits = c(min_lfc, max_lfc)) + 
  scale_color_manual(values = c("Differentially Expressed" = "#E64B35", 
                                "NS" = "#636363"),
                     labels = c("Differentially Expressed genes",
                                "Not significant")) 

coloredLung_Input <- ggplot(datalung_Input, aes(x = lfc, y = padj)) + 
  geom_point(aes(color = factor(color)), size = 1, alpha = 0.8, na.rm = T) + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5), legend.justification = c(0,1),
        legend.position = "none", panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour="black"),plot.margin = unit(c(0.1,0.1,0.2,0.1),"in")) +
  xlab(expression(log[2]("LPS" / "control"))) + 
  ylab(expression(-log[10]("q-value"))) + 
  geom_hline(yintercept = 1., colour = "black", linetype="dashed") + 
  geom_vline(xintercept = -1.0, colour = "black", linetype="dashed") + 
  geom_vline(xintercept = 1.0, colour = "black", linetype="dashed") +
  scale_y_continuous(limits = c(0, max_padj)) + 
  scale_x_continuous(limits = c(min_lfc, max_lfc)) + 
  scale_color_manual(values = c("Differentially Expressed" = "#E64B35", 
                                "NS" = "#636363"),
                     labels = c("Differentially Expressed genes",
                                "Not significant")) 

coloredLung_Bead <- ggplot(datalung_Bead, aes(x = lfc, y = padj)) + 
  geom_point(aes(color = factor(color)), size = 1, alpha = 0.8, na.rm = T) + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5), legend.justification = c(0,1),
        legend.position = "none", panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour="black"),plot.margin = unit(c(0.1,0.1,0.2,0.1),"in")) + 
  xlab(expression(log[2]("LPS" / "control"))) + 
  ylab(expression(-log[10]("q-value"))) + 
  geom_hline(yintercept = 1., colour = "black", linetype="dashed") + 
  geom_vline(xintercept = -1.0, colour = "black", linetype="dashed") + 
  geom_vline(xintercept = 1.0, colour = "black", linetype="dashed") +
  scale_y_continuous(limits = c(0, max_padj)) +
  scale_x_continuous(limits = c(min_lfc, max_lfc)) + 
  scale_color_manual(values = c("Differentially Expressed" = "#E64B35", 
                                "NS" = "#636363"),
                     labels = c("Differentially Expressed genes",
                                "Not significant")) 

blank<- rectGrob(gp=gpar(fill="white", col="white"))

title1 <- grobTree(rectGrob(gp=gpar(fill="white", col="white", lwd=5)),
                   textGrob("Brain", gp=gpar(col="black", fontsize= 13)))

title2 <- grobTree(rectGrob(gp=gpar(fill="white", col="white", lwd=5)),
                   textGrob("Heart", gp=gpar(col="black", fontsize= 13)))

title3 <- grobTree(rectGrob(gp=gpar(fill="white", col="white", lwd=5)),
                   textGrob("Kidney", gp=gpar(col="black", fontsize= 13)))

title4 <- grobTree(rectGrob(gp=gpar(fill="white", col="white", lwd=5)),
                   textGrob("Liver", gp=gpar(col="black", fontsize= 13)))

title5 <- grobTree(rectGrob(gp=gpar(fill="white", col="white", lwd=5)),
                   textGrob("Lung",  gp=gpar(col="black", fontsize= 13)))

title6 <- grobTree(rectGrob(gp=gpar(fill="white", col="white", lwd=5)),
                   textGrob("Tissue \nTranscriptome", rot = 90, gp=gpar(col="black", fontsize= 13)))

title7 <- grobTree(rectGrob(gp=gpar(fill="white", col="white", lwd=5)),
                   textGrob("EC \nTranslatome", rot = 90, gp=gpar(col="black", fontsize= 13)))

row1 <- ggarrange(blank, title1, title2, title3, title4, title5, ncol=6, widths = c(2, 12, 12, 12, 12, 12))
row2 <- ggarrange(title6, coloredBrain_Input, coloredHeart_Input, coloredKidney_Input, coloredLiver_Input, coloredLung_Input, 
                  ncol=6, nrow = 1, widths = c(2, 12, 12, 12, 12, 12))
row3 <- ggarrange(title7, coloredBrain_Bead, coloredHeart_Bead, coloredKidney_Bead, coloredLiver_Bead, coloredLung_Bead, 
                  ncol=6, nrow = 1, widths = c(2, 12, 12, 12, 12, 12))

figure <- ggarrange(row1,row2,row3,ncol=1, nrow=3, labels = c("","A", "B"),heights = c(1,10,10))


png("Figure7.png",width = 14, height = 6, units = "in", res=600)
figure
dev.off()

rm(list=ls()) 
