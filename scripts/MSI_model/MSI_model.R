# load libraries
library(caret)
library(doParallel)
library(stringr)
library(PRROC)
library(ROCit)
library(RColorBrewer)
library(org.Hs.eg.db)
source("model_functions.R")
set.seed(0202)

outputDir <- "../../results/MSI_model/"
dir.create(outputDir, recursive = T)

# read data
rna_data <- read.csv("../../results/RNA_Seq/genes/Counts/normalized_RNA_Seq.csv", header = T, row.names = 1)
genomic_metrics <- read.csv("../../results/sample_genome_metrics_matched.csv", header = T, row.names = 1)

# RNA gene annotation
gene_ids <- data.frame(ensembl = colnames(rna_data))

gene_ids$entrez <- mapIds(org.Hs.eg.db,
                          keys=as.character(gene_ids$ensembl), 
                          column="ENTREZID",
                          keytype="ENSEMBL",
                          multiVals="first")

gene_ids$symbol <- mapIds(org.Hs.eg.db,
                          keys=as.character(gene_ids$ensembl), 
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")

gene_ids$gene_description <- mapIds(org.Hs.eg.db,
                                    keys=as.character(gene_ids$ensembl), 
                                    column="GENENAME",
                                    keytype="ENSEMBL",
                                    multiVals="first")

# filter out genes with no Entrez ID
gene_ids <- gene_ids[!is.na(gene_ids$entrez),]

rna_data <- rna_data[,gene_ids$ensembl]

# annotate MSI and POLE/D1 as one variable
genomic_metrics$msi_pole <- as.character(genomic_metrics$msi_status)
genomic_metrics[which(genomic_metrics$msi_pole %in% c("MSI-L", "Indeterminate")), "msi_pole"] <- "MSS"
genomic_metrics[is.na(genomic_metrics$msi_pole), "msi_pole"] <- "-"
genomic_metrics[which((genomic_metrics$POLE_mut == T | genomic_metrics$POLD1_mut == T) & genomic_metrics$msi_pole != "MSI-H"), "msi_pole"] <- "POLE/D1"
genomic_metrics$msi_pole <- as.factor(genomic_metrics$msi_pole)

# remove duplicated sample ids
genomic_filtered <- genomic_metrics[!duplicated(genomic_metrics$sample_ids),]

rownames(genomic_filtered) <- genomic_filtered$sample_ids

rna_data <- rna_data[rownames(genomic_filtered),]

genomic_filtered$msi_state <- "MSS"
genomic_filtered[genomic_filtered$msi_status %in% c("MSI-H"), "msi_state"] <- "MSI"
genomic_filtered$msi_state <- factor(genomic_filtered$msi_state)

# create labels for model
gene_bp_status <- factor(genomic_filtered[match(rownames(rna_data), genomic_filtered$sample_ids),]$msi_state)
names(gene_bp_status) <- rownames(rna_data)

set.seed(203)

# run model training only if there are enough bps
splits <- createDataPartition(gene_bp_status, times = 1, p = 0.65, list = TRUE)
train_labels <- gene_bp_status[splits$Resample1]
test_labels <- gene_bp_status[-splits$Resample1]

# subsample rna-seq data
RNA_train <- rna_data[names(train_labels),]
RNA_test <- rna_data[names(test_labels),]
  
results <- run_model(RNA_train, train_labels, method = "rf", folds = 5,
                     tune_length = 10, features = 100, 
                     seq_length = 25, iteration_length = 5,
                     cores = 30, min_features = 1, 
                     major_class = "MSI", minor_class = "MSS")

ggplot(results$best_features) + xlim(10, 100)
ggplot(results$best_model)

prediction <- predict(results$best_model, RNA_test)
conf_result <- confusionMatrix(data = prediction, reference = test_labels)

prediction <- predict(results$best_model, RNA_test, type = "prob")

roc_bin <- rocit(score = prediction$MSI, class = test_labels, negref = "MSS")

balanced_test <- RNA_test
balanced_test$label <- test_labels
balanced_test <- sampling::strata(balanced_test, stratanames = "label", size = c(24,24))
balanced_test_data <- RNA_test[balanced_test$ID_unit,]

balanced_prediction <- predict(results$best_model, balanced_test_data)
conf_result <- confusionMatrix(data = balanced_prediction, reference = balanced_test$label)

balanced_prediction <- predict(results$best_model, balanced_test_data, type = "prob") 


plot(balanced_roc, YIndex = F, 
     values = F, legend = F)
legend("bottomright", col = c(1,alpha("grey50")),
       c(paste("Empirical ROC curve (AUC: ", format(balanced_roc$AUC, digits = 3), ")", sep=""), "Chance line"), 
       lty= c("solid", "dashed"), lwd = 2)

library(pROC)
myRoc <- roc(predictor = balanced_prediction$MSI, response = balanced_test$label, positive = "MSI")
plot(myRoc, print.auc = T, print.auc.x = 0.2, print.auc.y = 0.2, las=1, bty="l",
     cex.lab = 1.5, xlim = c(1,0), ylim = c(0,1))

roc_data <- data.frame(obs = balanced_prediction, class = balanced_test$label)
roc_data$class <- factor(roc_data$class)

library(plotROC)
pdf(paste(outputDir, "MSI", "_ROC.pdf", sep=""), width = 7, height = 7)
g <- ggplot(roc_data, aes(m=obs.MSS, d=factor(class, levels = c("MSI", "MSS")))) +
  geom_roc(n.cuts=0) +
  coord_equal() + 
  style_roc(xlab = "1 - Specificity (FPR)", ylab = "Sensitivity (TPR)", 
            theme = theme_classic, major.breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  geom_abline(linetype = "dashed") + theme_classic(base_size = 21) + 
  theme(axis.text = element_text(size = 18))

g + annotate("text", x=0.75, y=0.25, 
             label=paste("AUC =", round((calc_auc(g))$AUC, 2)),
             size = 6)
dev.off()

save.image(paste(outputDir, "MSI_model.RData", sep=""))
