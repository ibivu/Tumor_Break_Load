source("../Break_calling/cnd_functions.R")
library(dplyr)

# read filtered breakpoint calls
outputdir <- "../../results/SNP6/masked_copy_number_segment/primary_tumor/"
filtered_bp_Calls <- read.csv(paste(outputdir, "filterd_tumor_breakpoint_data.csv", sep=""), header = T)
all_bp_calls <- read.csv(paste(outputdir, "unfiltered_tumor_breakpoint_data.csv", sep=""), header = T)

# Count ammount of Somatic Breakpoints for every sample
Tumor_Break_Load <- bp_counting(filtered_bp_Calls, 
                                no_bp_sample = setdiff(all_bp_calls$sample, filtered_bp_Calls$sample))

# write TBL
write.csv(Tumor_Break_Load, paste(outputdir, "tumor_break_load.csv", sep = ""))