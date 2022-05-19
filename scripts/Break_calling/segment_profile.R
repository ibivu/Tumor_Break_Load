# create example copy number profiles (figure 1)
source("cnd_functions.R")
library(karyoploteR)
library(CopyNumberPlots)
library(scales)
library(circlize)


# read copy number data (set data folder)
tumor_workdir <- '../../data/SNP6/masked_copy_number_segment/primary_tumor/data'

# load segments
segments <- load_segment_files(tumor_workdir)

gg <- seqlengths(filterChromosomes(getGenome("hg38")))

# create segment grange
grange_segments <- GRanges(seqnames = segments$Chromosome, 
                           ranges = IRanges(start = segments$Start, end = segments$End), 
                           y = segments$Segment_Mean, sample = segments$GDC_Aliquot)

# change grange style to match chromosome name style
seqlevelsStyle(grange_segments) <- "UCSC"

# copy number profile plot function
copy_profile <- function(grange_segments, sample_id) {
  grange_sample <- grange_segments[grange_segments$sample == sample_id]
  
  grange_sample$lrr <- grange_sample$y
  grange_sample$cn <- grange_sample$y
  
  # cap plotting to same max values to have plots consistent
  ymax = 2
  ymin = -2 
  
  pp <- getDefaultPlotParams(4)
  pp$leftmargin = 0.1
  pp$rightmargin = 0.01
  pp$ideogramheight = 5
  kp <- plotKaryotype(genome = "hg38", plot.type = 4, cex = 1.5, 
                      labels.plotter = NULL,
                      chromosomes = as.character(unique(grange_sample@seqnames)), 
                      plot.params = pp)

  evens <- ((1:length(kp$chromosomes))%%2)==0
  
  chr.names <- kp$chromosomes
  even.names <- chr.names
  even.names[!evens] <- ""
  
  odd.names <- chr.names
  odd.names[evens] <- ""
  
  kpAddChromosomeNames(kp, chr.names = odd.names, yoffset = 10, cex=1.2)
  kpAddChromosomeNames(kp, chr.names = even.names, yoffset = -1, cex = 1.2)
  kpAbline(kp, h=0, ymin=ymin, ymax=ymax, lty=2, col="grey")
  pos_grange <- grange_sample[grange_sample$y > 0]
  neg_grange <- grange_sample[grange_sample$y < 0]

  plotCopyNumberCallsAsLines(kp, cn.calls = grange_sample, 
                             ymin = ymin, ymax = ymax, lwd=4, labels = NA, 
                             add.axis = F, style = "segments")
  kpAddChromosomeSeparators(kp, lwd=2, col="#666666")
  kpAxis(kp, ymin = ymin, ymax = ymax, tick.pos = seq(ymin, ymax, 1), cex=2.2)
  kpAddLabels(kp, labels = expression('Log'[2]*' ratio'), cex = 2.5, srt=90, 
              pos=3, label.margin = 0.065)
}

dir.path = "../../results/segment_profile/"
dir.create(dir.path)

samples <- unique(grange_segments$sample)
# create copy number plot for every sample
for (sample in samples) {
  pdf(paste(dir.path, "copy_number_profile_", sample, ".pdf", sep=""), width = 15, height = 7)
  copy_profile(grange_segments, sample)
  dev.off()
}

pdf(paste(dir.path, "copy_number_profile_", "c33148e3-a93e-4631-a034-d5302409ad65", "_whole_genome.pdf", sep=""), width = 15, height = 7)

# create example profiles of the SAS and BS with a genome wide view and a zoom in
grange_sample <- grange_segments[grange_segments$sample == "c33148e3-a93e-4631-a034-d5302409ad65",]
grange_sample

grange_sample$lrr <- grange_sample$y
grange_sample$cn <- grange_sample$y

ymax = 1
ymin = -1 

pp <- getDefaultPlotParams(4)
pp$leftmargin = 0.1
pp$rightmargin = 0.01
pp$ideogramheight = 5
kp <- plotKaryotype(genome = "hg38", plot.type = 4, cex = 1.5, 
                    labels.plotter = NULL,
                    chromosomes = as.character(unique(grange_sample@seqnames)), 
                    plot.params = pp)
evens <- ((1:length(kp$chromosomes))%%2)==0

chr.names <- kp$chromosomes
even.names <- chr.names
even.names[!evens] <- ""

odd.names <- chr.names
odd.names[evens] <- ""

kpAddChromosomeNames(kp, chr.names = odd.names, yoffset = 10, cex=1.2)
kpAddChromosomeNames(kp, chr.names = even.names, yoffset = -1, cex = 1.2)
kpAbline(kp, h=0, ymin=ymin, ymax=ymax, lty=2, col="grey")
pos_grange <- grange_sample[grange_sample$y > 0]
neg_grange <- grange_sample[grange_sample$y < 0]

plotCopyNumberCallsAsLines(kp, cn.calls = grange_sample, 
                           ymin = ymin, ymax = ymax, lwd=4, labels = NA, 
                           add.axis = F, style = "segments")
kpAddChromosomeSeparators(kp, lwd=2, col="#666666")
kpAxis(kp, ymin = ymin, ymax = ymax, tick.pos = seq(ymin, ymax, 1), cex=2.2)
kpAddLabels(kp, labels = expression('Log'[2]*' ratio'), cex = 2.5, srt=90, 
            pos=3, label.margin = 0.065)
dev.off()

# zoom in on chromosome 1
pdf(paste(dir.path, "copy_number_profile_", "c33148e3-a93e-4631-a034-d5302409ad65", "_zoom_in.pdf", sep=""), width = 10, height = 7)
grange_sample <- grange_segments[grange_segments$sample == "c33148e3-a93e-4631-a034-d5302409ad65" & seqnames(grange_segments) == "chr1",]
grange_sample

grange_sample$lrr <- grange_sample$y
grange_sample$cn <- grange_sample$y

ymax = 1
ymin = -1 

pp <- getDefaultPlotParams(4)
pp$leftmargin = 0.1
pp$rightmargin = 0.01
pp$ideogramheight = 5
kp <- plotKaryotype(genome = "hg38", plot.type = 4, cex = 1.5, 
                    labels.plotter = NULL,
                    chromosomes = as.character(unique(grange_sample@seqnames)), 
                    plot.params = pp)
#kpAddChromosomeNames(kp, srt=45, cex=2.2)#, yoffset = 5)
evens <- ((1:length(kp$chromosomes))%%2)==0

chr.names <- kp$chromosomes
even.names <- chr.names
even.names[!evens] <- ""

odd.names <- chr.names
odd.names[evens] <- ""

kpAddChromosomeNames(kp, chr.names = odd.names, yoffset = 10, cex=1.2)
kpAddChromosomeNames(kp, chr.names = even.names, yoffset = -1, cex = 1.2)
kpAbline(kp, h=0, ymin=ymin, ymax=ymax, lty=2, col="grey")
pos_grange <- grange_sample[grange_sample$y > 0]
neg_grange <- grange_sample[grange_sample$y < 0]

plotCopyNumberCallsAsLines(kp, cn.calls = grange_sample, 
                           ymin = ymin, ymax = ymax, lwd=4, labels = NA, 
                           add.axis = F, style = "segments")
kpAddChromosomeSeparators(kp, lwd=2, col="#666666")
kpAxis(kp, ymin = ymin, ymax = ymax, tick.pos = seq(ymin, ymax, 1), cex=2.2)
kpAddLabels(kp, labels = expression('Log'[2]*' ratio'), cex = 2.5, srt=90, 
            pos=3, label.margin = 0.065)
dev.off()


pdf(paste(dir.path, "copy_number_profile_", "153dfce3-3397-4e87-ba14-4c69c9c909df", "_whole_genome.pdf", sep=""), width = 12, height = 7)
# create example profiles of the SAS and BS with a genome wide view and a zoom in
grange_sample <- grange_segments[grange_segments$sample == "153dfce3-3397-4e87-ba14-4c69c9c909df",]
grange_sample

grange_sample$lrr <- grange_sample$y
grange_sample$cn <- grange_sample$y

ymax = 1
ymin = -1 

pp <- getDefaultPlotParams(4)
pp$leftmargin = 0.1
pp$rightmargin = 0.01
pp$ideogramheight = 5
kp <- plotKaryotype(genome = "hg38", plot.type = 4, cex = 1.5, 
                    labels.plotter = NULL,
                    chromosomes = as.character(unique(grange_sample@seqnames)), 
                    plot.params = pp)
#kpAddChromosomeNames(kp, srt=45, cex=2.2)#, yoffset = 5)
evens <- ((1:length(kp$chromosomes))%%2)==0

chr.names <- kp$chromosomes
even.names <- chr.names
even.names[!evens] <- ""

odd.names <- chr.names
odd.names[evens] <- ""

kpAddChromosomeNames(kp, chr.names = odd.names, yoffset = 10, cex=1.2)
kpAddChromosomeNames(kp, chr.names = even.names, yoffset = -1, cex = 1.2)
kpAbline(kp, h=0, ymin=ymin, ymax=ymax, lty=2, col="grey")
pos_grange <- grange_sample[grange_sample$y > 0]
neg_grange <- grange_sample[grange_sample$y < 0]

plotCopyNumberCallsAsLines(kp, cn.calls = grange_sample, 
                           ymin = ymin, ymax = ymax, lwd=4, labels = NA, 
                           add.axis = F, style = "segments")
kpAddChromosomeSeparators(kp, lwd=2, col="#666666")
kpAxis(kp, ymin = ymin, ymax = ymax, tick.pos = seq(ymin, ymax, 1), cex=2.2)
kpAddLabels(kp, labels = expression('Log'[2]*' ratio'), cex = 2.5, srt=90, 
            pos=3, label.margin = 0.065)
dev.off()

# zoom in on chromosome 13
pdf(paste(dir.path, "copy_number_profile_", "153dfce3-3397-4e87-ba14-4c69c9c909df", "_zoom_in.pdf", sep=""), width = 12, height = 7)
grange_sample <- grange_segments[grange_segments$sample == "153dfce3-3397-4e87-ba14-4c69c9c909df" & seqnames(grange_segments) == "chr13",]
grange_sample

grange_sample$lrr <- grange_sample$y
grange_sample$cn <- grange_sample$y

ymax = 1
ymin = -1 

pp <- getDefaultPlotParams(4)
pp$leftmargin = 0.1
pp$rightmargin = 0.01
pp$ideogramheight = 5
kp <- plotKaryotype(genome = "hg38", plot.type = 4, cex = 1.5, 
                    labels.plotter = NULL,
                    chromosomes = as.character(unique(grange_sample@seqnames)), 
                    plot.params = pp)
#kpAddChromosomeNames(kp, srt=45, cex=2.2)#, yoffset = 5)
evens <- ((1:length(kp$chromosomes))%%2)==0

chr.names <- kp$chromosomes
even.names <- chr.names
even.names[!evens] <- ""

odd.names <- chr.names
odd.names[evens] <- ""

kpAddChromosomeNames(kp, chr.names = odd.names, yoffset = 10, cex=1.2)
kpAddChromosomeNames(kp, chr.names = even.names, yoffset = -1, cex = 1.2)
kpAbline(kp, h=0, ymin=ymin, ymax=ymax, lty=2, col="grey")
pos_grange <- grange_sample[grange_sample$y > 0]
neg_grange <- grange_sample[grange_sample$y < 0]

plotCopyNumberCallsAsLines(kp, cn.calls = grange_sample, 
                           ymin = ymin, ymax = ymax, lwd=4, labels = NA, 
                           add.axis = F, style = "segments")
kpAddChromosomeSeparators(kp, lwd=2, col="#666666")
kpAxis(kp, ymin = ymin, ymax = ymax, tick.pos = seq(ymin, ymax, 1), cex=2.2)
kpAddLabels(kp, labels = expression('Log'[2]*' ratio'), cex = 2.5, srt=90, 
            pos=3, label.margin = 0.065)
dev.off()

