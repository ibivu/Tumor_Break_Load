library(dplyr)

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
    segs[1:(nrow(segs)-1),]$ID,        # Sample
    segs[1:(nrow(segs)-1),]$Chromosome,    # Chromosome
    segs[2:nrow(segs),]$ID,            # Sample + 1
    segs[2:nrow(segs),]$Chromosome,        # Chromosome + 1
    segs[1:(nrow(segs)-1),]$End,           # End, breakpoints x start (end segment 1)
    segs[2:nrow(segs),]$Start,             # +1 Start, breakpoint x end (start segment 2)
    segs[1:(nrow(segs)-1),]$Mean,  # (copy number segment 1)
    segs[2:nrow(segs),]$Mean,      # (copy number segment 2)
    segs[1:(nrow(segs)-1),]$Num_Markers,
    segs[2:nrow(segs),]$Num_Markers)             
  
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

calculate_amount_of_breakpoints_per_filtering <- function(segs){
  segment_filtering_values <- c(0, 10^seq(0, 4, 0.05)) # why this range?
  break_filtering_values <- c(0, 10^seq(-2, 1, 0.05)) # why this range?
  filter_count_matrix <- matrix(NA, length(segment_filtering_values), length(break_filtering_values))
  
  for(i in 1:length(segment_filtering_values)){
    print(i)
    # Filter segments, calculate breakpoints
    filtered_segs <- segs[segs$Num_Probes > segment_filtering_values[i],]
    filtered_breakpoints <- calculate_breakpoints(filtered_segs)
    
    # Loop through break sizes and fill matrix
    for(j in 1:length(break_filtering_values)){
      # Calculate amount of breakpoints at filtering, fill matrix
      amount_of_filter_breakpoints <- sum(filtered_breakpoints$break_size > break_filtering_values[j])
      filter_count_matrix[i,j] <- amount_of_filter_breakpoints
    }
  }
  return(filter_count_matrix)
}

make_filtered_plot <- function(segs, type, filter, outputdir, index){
  
  filtered_segs <- segs[segs$Num_Probes >= filter,]
  bps <- calculate_breakpoints(filtered_segs)
  
  pdf <- data.frame(bps$smallest_adjacent_segment,
                    bps$break_size)
  colnames(pdf) <- c('sas','bs')
  
  plottitle <- paste(type,', Segments smaller than ', round(filter), ' filtered out', sep='')
  plotfilename <- paste(type, index, '_filter_',round(filter),'_heatmap.png', sep='')
  plotpath <- paste(outputdir, plotfilename, sep='/')
  
  make_bp_heatmap(pdf, plottitle, plotpath)
}

make_bp_heatmap <- function(pdf, plottitle, plotpath){
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r <- rf(32)
  
  plt <- ggplot(pdf, aes(sas,bs))
  
  bp_heatmap <- plt + stat_bin2d(bins=100) + scale_fill_gradientn(colours=r, trans="log") + 
    scale_x_continuous(trans='log10', limits = c(1, 100000)) +
    scale_y_continuous(trans='log10', limits = c(0.01, 10)) +
    xlab('Smallest adjacent segment (# probes)') +
    ylab('Break size (delta log2ratio)') +
    theme(text = element_text(size=20)) +
    ggtitle(plottitle)
  
  png(filename = plotpath, width = 1000, height = 750, units = "px") 
  plot(bp_heatmap)
  dev.off()
}

filter_plots <- function(seg_data, type, outputdir){
  segment_filtering_values <- 10^seq(0, 4, 0.1)
  
  for(index_val in 1:length(segment_filtering_values)){
    
    print(index_val)
    
    filter <- segment_filtering_values[index_val]
    make_filtered_plot(seg_data, type, filter, outputdir, index_val)
  }
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