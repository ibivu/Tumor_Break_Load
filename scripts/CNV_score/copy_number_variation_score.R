source("../Break_calling/cnd_functions.R")
library(dplyr)
library(stringr)
library(GenomicRanges)
library(biomaRt)
library(ggplot2)

# read masked copy number variation
outputdirTumor <- "../../results/SNP6/masked_copy_number_segment/primary_tumor/"

dir.create(outputdirTumor, recursive = T)

# Set the directories for which you want to load the segments
tumor_workdir <- '../../data/SNP6/masked_copy_number_segment/primary_tumor/data'

# Load segment files combined in one data frame
tumor_segs <- load_segment_files(tumor_workdir)

tumor_metadata <- read.csv("../../data/SNP6/masked_copy_number_segment/primary_tumor/metadata.csv", header = T)

# convert file identifiers to TCGA identifiers
tumor_segs$submitter_id <- tumor_metadata[match(tumor_segs$GDC_Aliquot, tumor_metadata$entity_id), "entity_submitter_id"]

# create range and absolute segment mean
tumor_segs <- tumor_segs %>% mutate(range = End - Start, abs_cn = abs(Segment_Mean))

# indicate segment aberration state
tumor_segs$aberration <- F
tumor_segs[tumor_segs$abs_cn > 0.2, "aberration"] <- T
tumor_segs$aberration <- as.numeric(tumor_segs$aberration)

# aberrated segment size
tumor_segs$abr_seg <- tumor_segs$range * tumor_segs$aberration

sample_abr <- tumor_segs %>% group_by(submitter_id) %>% summarise(range = sum(range), abr_seg = sum(abr_seg))
sample_abr$fraction_abr <- sample_abr$abr_seg / sample_abr$range

# READ TBL data
Tumor_Break_Load <- read.csv(paste(outputdirTumor, "tumor_break_load.csv", sep=""), header = T, row.names = 1)

Tumor_Break_Load$fraction_abr <- sample_abr[match(Tumor_Break_Load$sample, sample_abr$submitter_id),]$fraction_abr

write.csv(Tumor_Break_Load, paste(outputdirTumor, "tumor_break_load.csv", sep=""))

# chromosomal arm aneuploidy
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
chromosome_annotation <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "strand", "band"),
              filters = "chromosome_name",
              values = levels(as.factor(tumor_segs$Chromosome)),
              mart = ensembl)

chromosome_annotation$arm <- str_replace(chromosome_annotation$band,
                                         pattern = "[0-9]+(\\.[0-9]+)?",
                                         replacement = "")

chromosome_annotation <- chromosome_annotation %>% group_by(chromosome_name, arm) %>% summarise(start = min(start_position), end = max(end_position))

chromosome_ranges <- GRanges(chromosome_annotation$chromosome_name, IRanges(chromosome_annotation$start, chromosome_annotation$end), arm=chromosome_annotation$arm)

tumor_ranges <- GRanges(tumor_segs$Chromosome, IRanges(start = tumor_segs$Start, end = tumor_segs$End),
                        submitter_id = as.character(tumor_segs$submitter_id), chr_range = tumor_segs$range, abr_seg = tumor_segs$abr_seg)

tumor_merged <- mergeByOverlaps(chromosome_ranges, tumor_ranges)
tumor_merged <- data.frame(tumor_merged)

chr_abr <- tumor_merged %>% group_by(submitter_id, tumor_ranges.seqnames, arm) %>% summarise(range = sum(chr_range), abr_seg = sum(abr_seg))
chr_abr$fraction_abr <- chr_abr$abr_seg / chr_abr$range

write.csv(chr_abr, paste(outputdirTumor, "chromsomal_aberration_score.csv", sep = ""))

# threshold aberration is 90%
chr_abr$aneuploidy_state <- F
chr_abr[chr_abr$fraction_abr >= 0.9, "aneuploidy_state"] <- T
chr_abr$submitter_id <- as.factor(chr_abr$submitter_id)

chr_score <- chr_abr %>% group_by(submitter_id) %>% summarise(aneuploidy_score = sum(aneuploidy_state))

Tumor_Break_Load$aneuploidy_score <- chr_score[match(Tumor_Break_Load$sample, chr_score$submitter_id),]$aneuploidy_score

ggplot(Tumor_Break_Load, aes(x=count, y=aneuploidy_score, color = msi_status)) + geom_point()
plot(Tumor_Break_Load$fraction_abr, Tumor_Break_Load$aneuploidy_score)

write.csv(Tumor_Break_Load, paste(outputdirTumor, "tumor_break_load.csv", sep=""))
