# load libraries
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(survival)
library(survminer)
library(pROC)

# clinical and break data
data <- read_excel("../data/1471-2407-14-121-S1.xls", col_names = T,  skip=3)

table(data[data$`PT_Stage TNM` %in% c(2,3),"Last clinical Status"])

# set recursion free survival events. Censor EVOL and Other events as there is no description of it
data$RFS_events <- NA
data[data$`Last clinical Status` == "R", "RFS_events"] <- 1
data[data$`Last clinical Status` == "OK", "RFS_events"] <- 0

# survival analysis
data$OS <- as.numeric(data$OS)
data$RFS <- as.numeric(data$RFS)
data$Event <- factor(data$Event, levels = c("Alive", "Dead"))
data$RFS_events <- as.numeric(data$RFS_events)

# convert month to days
data$RFS <- data$RFS*30.4368499

surv_data <- data
sfit = survfit(Surv(as.numeric(RFS), as.numeric(RFS_events))~1, data = surv_data)
plt = ggsurvplot(sfit)
plot.new()
print(plt, newpage = F)

# define TBL
surv_data$TBL <- surv_data$nb.BP

# determine the optimal cutoff with pROC coords
pdf("../results/TBL_threshold_selection_pROC.pdf")
plot.roc(surv_data$RFS_events, surv_data$TBL, thresholds = "best",
         print.thres="best")
dev.off()
tbl_threshold <- coords(roc(surv_data$RFS_events, surv_data$TBL), "best", 
                        ret="threshold")

# asign tbl classes based on the threshold
surv_data$TBL_class = "low"
surv_data[surv_data$TBL > tbl_threshold$threshold,"TBL_class"] = "high"

write.csv(surv_data, "../results/orsetti_surv_tbl_data.csv")

#tbl.cat <- surv_categorize(tbl.cut)
#tbl.cat

fit <- survfit(Surv(RFS, RFS_events) ~ TBL_class, data = surv_data)
ggsurvplot(
  fit,
  risk.table = T,
  pval = T,
  conf.int = F
)

# survival analysis (RFS) for high and low TBL quantiles (all stages)
sfit = survfit(Surv(RFS, RFS_events) ~ TBL_class, data = surv_data)
cox_TBL <- coxph(Surv(RFS, RFS_events) ~ factor(TBL_class, levels = c("low", "high")), 
                 data = surv_data)
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit)
coldata <- surv_data[!is.na(surv_data$RFS_events),]$TBL_class
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
    x = 3000, y = .3,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", 
                  round(exp(cox_TBL$coefficients), 3)," (",  
                  round(s_tbl$conf.int[3], 3), " - ", 
                  round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 5
  ) + xlab("Time (days)") + ylab("Relapse Free Survival (RFS)")
pdf("../results/RFS_survival_stage_1_4_tbl_thres_proc.pdf")
ggsurv
dev.off()
ggsurv

# survival analysis (RFS) for high and low TBL quantiles (stage 2 & 3)
sfit = survfit(Surv(RFS, RFS_events) ~ TBL_class, data = surv_data[surv_data$`PT_Stage TNM` %in% c(2,3),])
cox_TBL <- coxph(Surv(RFS, RFS_events) ~ factor(TBL_class, levels = c("low", "high")), 
                 data = surv_data[surv_data$`PT_Stage TNM` %in% c(2,3),])
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit)
coldata <- surv_data[!is.na(surv_data$RFS_events) & surv_data$`PT_Stage TNM` %in% c(2,3),]$TBL_class
my_xlab <- paste(levels(factor(coldata)), " (N=", table(coldata), ")", sep="")
ggsurv <- ggsurvplot(sfit, conf.int = F, pval = F, risk.table = F, 
                     palette = c("dodgerblue2", "orchid2"),
                     title = "Orsetti et al. MSS CRC stage 2 & 3",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot + theme(legend.text = element_text(size=12), 
                                   legend.title = element_text(size=14),
                                   plot.title = element_text(hjust = 0.5)) + 
  scale_color_discrete(name = "TBL", labels = my_xlab) +
  ggplot2::annotate(
    "text",
    x = 3000, y = .3,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", 
                  round(exp(cox_TBL$coefficients), 3)," (",  
                  round(s_tbl$conf.int[3], 3), " - ", 
                  round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 5
  ) + xlab("Time (days)") + ylab("Relapse Free Survival (RFS)")
pdf("../results/RFS_survival_stage_2_3_tbl_thres_proc.pdf")
ggsurv
dev.off()
ggsurv

# survival analysis (RFS) for high and low TBL quantiles (stage 2 & 3)
sfit = survfit(Surv(RFS, RFS_events) ~ TBL_class, data = surv_data[surv_data$`PT_Adjuvant Chemo` == "no" & 
                                                                     surv_data$`PT_Stage TNM` %in% c(2, 3),])
cox_TBL <- coxph(Surv(RFS, RFS_events) ~ factor(TBL_class, levels = c("low", "high")), 
                 data = surv_data[surv_data$`PT_Adjuvant Chemo` == "no" & 
                                    surv_data$`PT_Stage TNM` %in% c(2, 3),])
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit)
coldata <- surv_data[!is.na(surv_data$RFS_events) & surv_data$`PT_Stage TNM` %in% c(2,3) & 
                       surv_data$`PT_Adjuvant Chemo` == "no",]$TBL_class
my_xlab <- paste(levels(factor(coldata)), " (N=", table(coldata), ")", sep="")
ggsurv <- ggsurvplot(sfit, conf.int = F, pval = F, risk.table = F, 
                     palette = c("dodgerblue2", "orchid2"),
                     title = "Orsetti et al. MSS CRC stage 2 & 3",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot + theme(legend.text = element_text(size=12), 
                                   legend.title = element_text(size=14),
                                   plot.title = element_text(hjust = 0.5)) + 
  scale_color_discrete(name = "TBL", labels = my_xlab) +
  ggplot2::annotate(
    "text",
    x = 3000, y = .3,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", 
                  round(exp(cox_TBL$coefficients), 3)," (",  
                  round(s_tbl$conf.int[3], 3), " - ", 
                  round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 5
  ) + xlab("Time (days)") + ylab("Relapse Free Survival (RFS)")
pdf("../results/RFS_survival_stage_2_3_tbl_thres_proc_untreated.pdf")
ggsurv
dev.off()
ggsurv
