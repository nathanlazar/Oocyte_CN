# nathan dot lazar at gmail dot com

# Plots the ratio of observed over expected read counts in each of the 
# supplied bins for the whole sample as well as individual chromosomes
# are generated and placed in a directory with the base name of the 
# file plus _plots

# Usage: Rscript plot_bin_counts.R <bin_count_file.txt> <out_dir>

# Example: Rscript plot_bin_counts.R \
#          MAPPED_MASKED/10-24-14-B1_S7_L001.q30.rmdup.bin_200_counts.txt 

# Will place plots in MAPPED_MASKED/10-24-14-B1_S7_L001_plots

library(ggplot2)  # Used for plotting
library(dplyr)    # General data frame manipulation
library(DNAcopy)  # Implements circular binary segmentation
library(grid)     # Used for plotting 

# Set outlier ratios to mean ----------------------------------------------

rm_outliers <- function(counts, SD.cut=4, SD.shrink=2) {
  # Set ratio and count outliers (more than SD.cut s.d.s above mean of 
  # non-zero ratios) to the mean plus SD.shrink std. devs. of other ratios
  mean <- mean(counts$ratio[counts$count!=0])
  SD <- sd(counts$ratio[counts$count!=0])
  ratio.cut <- mean + SD.cut * SD
  counts$ratio[counts$ratio > ratio.cut] <- mean + SD.shrink * SD
  counts
}

# Use circular binary segmentation to call copy number states -------------

get_seg <- function(counts, base.cn) {
  counts.cna <- CNA(genomdat=counts$ratio * base.cn, 
                    chrom=counts$chr, 
                    maploc=1:nrow(counts), 
                    presorted=T)

  # This smoothes outliers. If a bin's ratio is more than 3 std. deviations
  # from it's nearest neighbor it will be shrunk to be 2 s.d. from the median.
  # The s.d. is taken from the distribution of the nearest 20 bins
  # Trim trims the most outlying 5% of data
  counts.smooth <- smooth.CNA(counts.cna, 
                              smooth.region=10, 
                              outlier.SD.scale=3,
                              smooth.SD.scale=2,
                              trim=0.05)

  # Run the circular binary segmentation.
  # Alpha is the significance level to accept change points.
  # nperm is the number of permutations run.
  # min.width is the smallest number of bins for a copy number segment
  # Prune undoes splits when the increase in sum of squares for unsplit 
  # is less than 5%
  counts.seg <- segment(counts.smooth, weights=rep(1,nrow(counts)), 
                        alpha = 1e-6, min.width=5,
                        nperm=1e6, undo.splits='sdundo',
                        undo.SD=.25, verbose=1)

#  plot(counts.seg, ylim=c(0,2))
  
  counts.seg$output$round <- round(counts.seg$output$seg.mean)
  counts.seg
}

# Recalibrate normalization based on the segmentation counts --------------

recalibrate <- function(counts, counts.seg, bin.size=200, base.cn=2, verbose=T) {
  # Assign putative copy number to each window
  counts$cn <- unlist(mapply(rep, counts.seg$output$round, counts.seg$output$num.mark))
  
  if(verbose) {
    print('Putative copy number means before correction:')
    counts %>% tbl_df %>% group_by(cn) %>% 
      summarise(mean.ratio=mean(ratio)) %>% print
  }
  
  # Correct for the differing amount of DNA in aneuploid samples
  control.tot <- sum(counts$count!=0) * bin.size * base.cn
  samp.tot <- sum(counts$cn[counts$count!=0])*(bin.size/base.cn) * 2
  # If all bins are called as copy number 0, don't correct
  if(samp.tot == 0) samp.tot <- control.tot
  counts$ratio <- counts$ratio * (samp.tot/control.tot) * base.cn
  
  if(verbose) {
    print('Putative copy number means after correction:')
    counts %>% tbl_df %>% group_by(cn) %>% 
      summarise(mean.ratio=mean(ratio)) %>% print
  }
  counts
}

# Main script -------------------------------------------------------------

# Load data
args <- commandArgs(trailingOnly = TRUE)
count.file <- args[1]

if(grepl('200', count.file)) bin.size=200
if(grepl('1000', count.file)) bin.size=1000

counts <- read.table(count.file, header=T, stringsAsFactors=F)
name <- sub('.bin_200_counts.txt', '', count.file)
name <- sub('.sort', '', name)
out.dir <- paste0(name, '_plots')
split.name <- strsplit(name, split='/', fixed=T)[[1]]
name <- split.name[length(split.name)]

# Make output directory if necessary 
if (!file.exists(out.dir)){
  dir.create(file.path(out.dir))
}

# Set bins w/ less than 2 reads to zero
read.cut <- 2
counts$count[counts$count < read.cut] <- 0
counts$ratio[counts$count < read.cut] <- 0
counts$per[counts$count < read.cut] <- 0

# Adjust expected counts removing bins w/ zero counts
counts$expected <- counts$expected / mean(counts$count != 0)
counts$ratio <- counts$per/counts$expected

# Shrink outliers above the mean by more than 3 std. devs
# to the mean plus 2 std. devs.
counts <- rm_outliers(counts, SD.cut=3, SD.shrink=2)

# If most bins are empty we expect a base copy number of 1 
# otherwise base copy number is set to 2
if(mean(counts$count == 0) > .5) {
  base.cn <- 1
} else base.cn <- 2

# Get segmentation
counts.seg <- get_seg(counts, base.cn)

# Recalibrate counts
counts <- recalibrate(counts, counts.seg, bin.size, base.cn, verbose=T)

# Prep for plotting
bins.per.chrom <- counts %>% tbl_df %>% group_by(chr) %>% summarise(len=length(chr)) 
bins.per.chrom$mid <- cumsum(bins.per.chrom$len) - (bins.per.chrom$len/2)
counts$idx <- 1:nrow(counts)
counts$bin <- unlist(sapply(bins.per.chrom$len, function(x) 1:x))
counts$col <- 2
counts$col[counts$chr %in% unique(counts$chr)[seq(1,21,2)]] <- 1
counts$col <- as.factor(counts$col)
ylim <- c(0, min(max(counts$ratio), 6))

# Plot whole genome
png(file=paste0(out.dir, '/all.png'), width = 480*2)
p <- ggplot(counts, aes(x=idx, y=ratio)) +
  geom_point(aes(colour=col)) +
  theme(legend.position="none") +
  labs(y='Copy Number', x='Chromosome') +
  scale_x_continuous(breaks=cumsum(bins.per.chrom$len), 
                     minor_breaks=NULL, labels=rep("", nrow(bins.per.chrom))) +
  scale_y_continuous(limits=ylim, breaks=0:6,
                     minor_breaks=NULL) +
  scale_colour_manual(values=c('#4861B3', '#70B333')) +
  theme(plot.title=element_text(size=24, face="bold"),
        axis.title=element_text(size=20),
        axis.text.x=element_text(size=14),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_text(size=20)) +
  geom_text(data=bins.per.chrom, size=7, colour='#7F7F7F',
            aes(x=bins.per.chrom$mid, vjust=2.5,
                y=rep(0,21), label=gsub("chr|chr0", "", unique(counts$chr)))) +
  labs(title=name) +
  geom_segment(data=counts.seg$output, 
               aes(x=loc.start, xend=loc.end, y=round, yend=round), size=2)

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)
dev.off()

for(chr in unique(counts$chr)) {
  chr.counts <- counts[counts$chr==chr,]
  chr.counts$idx <- chr.counts$bin
  chr.counts.seg <- counts.seg
  chr.counts.seg$output <- counts.seg$output[counts.seg$output$chrom==chr,]
  chr.start <- chr.counts.seg$output$loc.start[1]
  chr.counts.seg$output$loc.start <- chr.counts.seg$output$loc.start - chr.start
  chr.counts.seg$output$loc.end <- chr.counts.seg$output$loc.end - chr.start
  chr.bins.per.chrom <- bins.per.chrom[bins.per.chrom$chr==chr,]
  chr.bins.per.chrom$mid <- chr.bins.per.chrom$len/2
  chr.num <- gsub("chr|chr0", "", chr)

  # Hack to make colors alternate
  if(unique(chr.counts$col)==1) {
    cols <- c('#4861B3', '#70B333')
  } else cols <- c('#70B333', '#4861B3')
  
  # Make column showing breaks every 10 Mb
  chr.counts$mb <- chr.counts$end %/% 10000000
  mb10.bins <- chr.counts %>% tbl_df %>% group_by(mb) %>% summarise(mb.bins=length(idx))
  
  png(file=paste0(out.dir, '/', chr, '.png'), width = 480*2)
  p <- ggplot(chr.counts, aes(x=idx, y=ratio)) +
    geom_point(aes(colour=col)) +
    theme(legend.position="none") +
    labs(y='Copy Number', x="10Mb intervals") +
    scale_x_continuous(breaks=c(0,cumsum(mb10.bins$mb.bins)[-nrow(mb10.bins)]), 
                       labels=mb10.bins$mb) +
    scale_y_continuous(limits=ylim, breaks=0:6,
                       minor_breaks=NULL) +
    scale_colour_manual(values=cols) +
    theme(plot.title=element_text(size=24, face="bold"),
          axis.title=element_text(size=20),
          axis.text.x=element_text(size=20),
          axis.ticks.x=element_blank(), 
          axis.text.y=element_text(size=20)) +
    labs(title=paste(name, 'chromosome', chr.num)) +
    geom_segment(data=chr.counts.seg$output, aes(x=loc.start, xend=loc.end, 
                                      y=round, yend=round), size=2)
  print(p)
  dev.off()
}
