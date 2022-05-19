# assess the TBL threshold using the coord function from pROC
# Libraries
library(survival)
library(survminer)
library(plyr)
library(TCGAbiolinks)
library(stringr)
library(data.table)
library(tidyverse)
library(ggplot2)
library(lmvar)
library(caret)
library(ggpmisc)
library(MASS)
library(readxl)
library(pROC)

# results output directories
outputDir = "../../results/survival_analysis/"
dir.create(outputDir)

# read data, clinical data from liu et.al 2018 cell and genomic metrics
clinical_data = read_excel("../../data/mmc1.xlsx", sheet = "TCGA-CDR")
genome_metrics = read.csv("../../results/sample_genome_metrics_tbl_class.csv", 
                          header = T, row.names = 1)
full_metrics = read.csv("../../results/sample_genome_metrics_tbl_class_all.csv", 
                        header = T, row.names = 1)

# process survival data
survival_data = subset(clinical_data[which(clinical_data$bcr_patient_barcode %in% unique(full_metrics$barcode)),], 
                       select = c(bcr_patient_barcode, OS.time, OS, DSS, DSS.time, 
                                  DFI, DFI.time, PFI, PFI.time))

survival_data = survival_data[!duplicated(survival_data$bcr_patient_barcode),]
survival_data$OS.time <- as.numeric(as.character(survival_data$OS.time))
survival_data$DSS.time <- as.numeric(as.character(survival_data$DSS.time))
survival_data$DFI.time <- as.numeric(as.character(survival_data$DFI.time))
survival_data$PFI.time <- as.numeric(as.character(survival_data$PFI.time))
survival_data[survival_data == "[Not Available]"] <- NA
survival_data$OS <- as.numeric(survival_data$OS)
survival_data$DSS <- as.numeric(survival_data$DSS)
survival_data$DFI <- as.numeric(survival_data$DFI)
survival_data$PFI <- as.numeric(survival_data$PFI)

# overall survival plot of the TCGA
sfit = survfit(Surv(OS.time, OS)~1, data = survival_data)
plt = ggsurvplot(sfit)
plot.new()
print(plt, newpage = F)

# add TBL to survival data
survival_data$TBL = full_metrics[match(survival_data$bcr_patient_barcode, full_metrics$barcode), "count"]
survival_filtered = survival_data[match(genome_metrics$barcode, survival_data$bcr_patient_barcode),]

# determine the optimal cutoff with pROC coords
pdf(paste(outputDir, "TBL_threshold_selection_pROC.pdf", sep=""))
plot.roc(survival_filtered$DFI, survival_filtered$TBL, thresholds = "best",
         print.thres="best")
dev.off()
tbl_threshold <- coords(roc(survival_filtered$DFI, survival_filtered$TBL), "best", 
       ret="threshold")

# asign tbl classes based on the threshold
survival_filtered$TBL_class = "low"
survival_filtered[survival_filtered$TBL > tbl_threshold$threshold,"TBL_class"] = "high"

# retrieve TBL class labels from the threshold
tbl.cat <- survival_filtered[c("DFI.time", "DFI", "TBL", "TBL_class", "bcr_patient_barcode")]
tbl.cat

write.csv(tbl.cat, paste(outputDir, "dichotomized_TBL_pROC.csv", sep=""))

fit <- survfit(Surv(DFI.time, DFI) ~ TBL_class, data = tbl.cat)
ggsurvplot(
  fit,
  risk.table = T,
  pval = T,
  conf.int = F
)

# annotate with staging
query = GDCquery(project = c("TCGA-COAD", "TCGA-READ"),
                 data.category = "Clinical",
                 file.type = "xml",
                 barcode = full_metrics$barcode)
GDCdownload(query)
clinical = GDCprepare_clinic(query, clinical.info = "stage_event")

full_metrics$pathological_stage = clinical[match(full_metrics$barcode, 
                                                 clinical$bcr_patient_barcode), 
                                           "pathologic_stage"]

full_metrics$pathological_stage = factor(str_replace(full_metrics$pathological_stage, 
                                                     pattern = "[ABC]", replacement = ""))

full_filtered = full_metrics[full_metrics$pathological_stage %in% 
                               c("Stage I", "Stage II", "Stage III", "Stage IV"),]

full_filtered = full_filtered[!duplicated(full_filtered$sample_ids),]

rownames(full_filtered) = full_filtered$sample_ids

survival_filtered = survival_data[survival_data$bcr_patient_barcode %in% 
                                    full_filtered$barcode,]
survival_filtered$stage <- full_filtered[match(survival_filtered$bcr_patient_barcode, 
                                               full_filtered$barcode), 
                                         "pathological_stage"]

tbl.cat$stage <- survival_filtered[match(tbl.cat$bcr_patient_barcode, 
                                         survival_filtered$bcr_patient_barcode),]$stage

# survival analysis (DFS) for high and low TBL quantiles (all stages)
sfit = survfit(Surv(DFI.time, DFI) ~ TBL_class, data = tbl.cat)
cox_TBL <- coxph(Surv(DFI.time, DFI) ~ factor(TBL_class, levels = c("low", "high")), data = tbl.cat)
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit)
coldata <- tbl.cat[!is.na(tbl.cat$DFI),]$TBL_class
my_xlab <- paste(levels(factor(coldata)), " (N=", table(coldata), ")", sep="")
ggsurv <- ggsurvplot(sfit, conf.int = F, pval = F, risk.table = F, 
                     palette = c("dodgerblue2", "orchid2"),
                     title = "Kaplan-Meier Curve for TBL colorectal cancer Stage 1-4",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot + theme(legend.text = element_text(size=12), legend.title = element_text(size=14),
                                   plot.title = element_text(hjust = 0.5)) + 
  scale_color_discrete(name = "TBL", labels = my_xlab) +
  ggplot2::annotate(
    "text",
    x = 4000, y = .3,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", round(exp(cox_TBL$coefficients), 3)," (",  round(s_tbl$conf.int[3], 3), " - ", round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 5
  ) + xlab("Time (Days)") + ylab("Disease Free Interval (DFI)")
pdf(paste(outputDir, "survival_analysis_DFI_MSS_all_thres_proc.pdf", sep=""))
ggsurv
dev.off()
ggsurv

sfit = survfit(Surv(DFI.time, DFI) ~ TBL_class, data = tbl.cat[tbl.cat$stage %in% c("Stage II", "Stage III"),])
cox_TBL <- coxph(Surv(DFI.time, DFI) ~ factor(TBL_class, levels = c("low", "high")), 
                 data = tbl.cat[tbl.cat$stage %in% c("Stage II", "Stage III"),])
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit)
coldata <- tbl.cat[!is.na(tbl.cat$DFI) & tbl.cat$stage %in% c("Stage II", "Stage III"),]$TBL_class
my_xlab <- paste(levels(factor(coldata)), " (N=", table(coldata), ")", sep="")
ggsurv <- ggsurvplot(sfit, conf.int = F, pval = F, risk.table = F, 
                     palette = c("dodgerblue2", "orchid2"),
                     title = "TCGA MSS CRC Stage 2 & 3",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot + theme(legend.text = element_text(size=12), legend.title = element_text(size=14),
                                   plot.title = element_text(hjust = 0.5)) + 
  scale_color_discrete(name = "TBL", labels = my_xlab) +
  ggplot2::annotate(
    "text",
    x = 4000, y = .3,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", round(exp(cox_TBL$coefficients), 3)," (",  
                  round(s_tbl$conf.int[3], 3), " - ", round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 5
  ) + xlab("Time (Days)") + ylab("Disease Free Survival (DFS)")
pdf(paste(outputDir, "survival_analysis_DFI_MSS_23_thres_proc.pdf", sep=""))
ggsurv
dev.off()
ggsurv
