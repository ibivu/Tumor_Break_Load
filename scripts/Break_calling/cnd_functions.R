# Function to load segment files from within a folder or nested folder structure and combine them into one data frame
load_segment_files <- function(workdir, pattern = F){
  # extract segment files recursively from folders in workdir 
  seg_files <- list.files(workdir, 
                          full.names = T, recursive = T)[grep('seg.v2.txt$',list.files(workdir, recursive = T))]
  all_segs <- data.frame()
  
  # merge segment file dataframes
  for(i in 1:length(seg_files)){
    #print(i)
    if (pattern != F) {
      if (!grepl(pattern, seg_files[i])) {
        segs_of_index <- read.table(seg_files[i], head=T)
        segs_of_index$filename <- basename(seg_files[i])
        all_segs <- rbind(all_segs, segs_of_index)
      }
    } else {
      segs_of_index <- read.table(seg_files[i], head=T)
      segs_of_index$filename <- basename(seg_files[i])
      all_segs <- rbind(all_segs, segs_of_index)
    }
  }
  return(all_segs)
}

# Calculate breakpoint metrics from a (filtered) segments data frame.
calculate_breakpoints <- function(segs){
  
  bp_construction_df <- data.frame(
    segs[1:(nrow(segs)-1),]$submitter_id,        # Sample
    segs[1:(nrow(segs)-1),]$Chromosome,    # Chromosome
    segs[2:nrow(segs),]$submitter_id,            # Sample + 1
    segs[2:nrow(segs),]$Chromosome,        # Chromosome + 1
    segs[1:(nrow(segs)-1),]$End,           # End, breakpoints x start (end segment 1)
    segs[2:nrow(segs),]$Start,             # +1 Start, breakpoint x end (start segment 2)
    segs[1:(nrow(segs)-1),]$Segment_Mean,  # (copy number segment 1)
    segs[2:nrow(segs),]$Segment_Mean,      # (copy number segment 2)
    segs[1:(nrow(segs)-1),]$Num_Probes,
    segs[2:nrow(segs),]$Num_Probes)             
  
  colnames(bp_construction_df) <- c('sample','chromosome','sample_1','chromosome_1','x_start','x_end','y_start','y_end')
  
  # remove overlap between chromosomes or files
  bp_df <- bp_construction_df[bp_construction_df$sample == bp_construction_df$sample_1 & 
                                bp_construction_df$chromosome == bp_construction_df$chromosome_1,] 
  
  bp_df <- bp_df[,-c(3,4)] # Remove the +1 rows
  
  # Calculate new info (e.g. x and y break sizes)
  bp_df$dx <- bp_df$x_end - bp_df$x_start # size between the segments
  bp_df$xpos <- (bp_df$x_end + bp_df$x_start)/2 # middle pos of the break
  bp_df$dy <- bp_df$y_end - bp_df$y_start # Segment_mean difference
  bp_df$break_size <- abs(bp_df$dy) # break size
  
  bp_df$smallest_adjacent_segment <- apply(bp_df[,7:8],1,min) # select lowest number of probes for the adjacent segements
  bp_df <- bp_df[,-c(7,8)] # Remove the +1 rows
  
  return(bp_df)
}

bp_counting <- function(data, sample_column="sample", no_bp_sample=NULL){
  # Count chromosomal breaks per sample
  # data: breakpoint call dataframe from calculate_breakpoints function
  # sample_column: string of column name with sample identifiers
  # no_bp_sample: vector of sample identifiers not included in the data with 0 breaks
  count_data <- data %>% group_by(.dots = as.symbol(sample_column)) %>% summarise(count = n())
  if (length(no_bp_sample) != 0) {
    m <- matrix(0, ncol = dim(count_data)[2], nrow = length(no_bp_sample))
    m <- as.data.frame(m)
    colnames(m) <- colnames(count_data)
    m$sample <- no_bp_sample
    
    count_data <- rbind(count_data, m)
  }
  return(count_data)
}

barcode <- function(x){
  # Convert TCGA identifier to TCGA barcode (TCGA case identifier)
  return(paste(strsplit(as.character(x), "-")[[1]][1:3], collapse="-"))
}