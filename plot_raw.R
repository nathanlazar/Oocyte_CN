# nathan dot lazar at gmail dot com

# Plots the read counts in each of the supplied bins for the
# whole sample as well as individual chromosomes. png files
# are generated and placed in a directory with the base name of the 
# file plus _plots

# Usage: Rscript plot_raw.R <bin_count_file.cn> <out_dir>

# Example: Rscript plot_raw.R COUNTS/4000/rh150409-1-b1-B1_S1.cn'

library(dplyr)    # data frame manipulation
library(ggplot2)  # Used for plotting
library(grid)     # Used for plotting 

# Main script -------------------------------------------------------------

# Load data
args <- commandArgs(trailingOnly = TRUE)
cn.file <- args[1]
out.dir <- sub('.cn', '_plots', cn.file)

name <- sub('.cn', '', cn.file)
split.name <- strsplit(name, split='/', fixed=T)[[1]]
name <- split.name[length(split.name)]

counts <- read.table(cn.file, header=T, stringsAsFactors=F)
                     
# Make output directory if necessary 
if (!file.exists(out.dir))
  dir.create(file.path(out.dir))

# Make a dataframe of the copy number calls
# counts.seg <- tbl_df(counts) %>% 
#   select(chr, idx, cn) %>%
#   group_by(chr, cn) %>%
#   summarise(start=min(idx), end=max(idx))


counts$chr <- factor(counts$chr, levels=c(paste0('chr', 1:20), 'chrX'))

# Prep for plotting
bins.per.chrom <- counts %>% tbl_df %>% group_by(chr) %>% summarise(len=length(chr)) 
bins.per.chrom$mid <- cumsum(bins.per.chrom$len) - (bins.per.chrom$len/2)
n.chr <- length(unique(counts$chr))
counts$col <- 2
counts$col[counts$chr %in% unique(counts$chr)[seq(1,n.chr,2)]] <- 1
counts$col <- as.factor(counts$col)
ylim <- c(0, max(counts$count))

# Plot whole genome
png(file=paste0(out.dir, '/all_raw.png'), width = 480*2)
p <- ggplot(counts, aes(x=idx, y=count)) +
  geom_point(aes(colour=col)) +
  theme(legend.position="none") +
  labs(y='Reads in each bin', x='Chromosome') +
  scale_x_continuous(breaks=cumsum(bins.per.chrom$len), 
                     minor_breaks=NULL, labels=rep("", nrow(bins.per.chrom))) +
  scale_y_continuous(limits=ylim, minor_breaks=NULL) +
  scale_colour_manual(values=c('#4861B3', '#70B333')) +
  theme(plot.title=element_text(size=24, face="bold", hjust=.5),
        axis.title=element_text(size=20),
        axis.text.x=element_text(size=14),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_text(size=20)) +
  geom_text(data=bins.per.chrom, size=6, colour='#7F7F7F',
            aes(x=bins.per.chrom$mid, vjust=2.5,
                y=rep(0,nrow(bins.per.chrom)), 
                label=sub('chr', '', bins.per.chrom$chr))) +
  labs(title=name) 


gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)
dev.off()

for(chr in unique(counts$chr)) {
  chr.counts <- counts[counts$chr==chr,]
  chr.counts$idx <- chr.counts$bin
  chr.bins.per.chrom <- bins.per.chrom[bins.per.chrom$chr==chr,]
  chr.bins.per.chrom$mid <- chr.bins.per.chrom$len/2
  chr.num <- gsub("chr", "", chr)

  # Hack to make colors alternate
  if(unique(chr.counts$col)==1) {
    cols <- c('#4861B3', '#70B333')
  } else cols <- c('#70B333', '#4861B3')
  
  # Make column showing breaks every 10 Mb
  chr.counts$mb <- chr.counts$end %/% 10000000
  mb10.bins <- chr.counts %>% tbl_df %>% group_by(mb) %>% summarise(mb.bins=length(idx))
  
  png(file=paste0(out.dir, '/', chr, '_raw.png'), width = 480*2)
  p <- ggplot(chr.counts, aes(x=idx, y=count)) +
    geom_point(aes(colour=col)) +
    theme(legend.position="none") +
    labs(y='Reads in each bin', x="10Mb intervals") +
    scale_x_continuous(breaks=c(0,cumsum(mb10.bins$mb.bins)[-nrow(mb10.bins)]), 
                       labels=mb10.bins$mb) +
    scale_y_continuous(limits=ylim, minor_breaks=NULL) +
    scale_colour_manual(values=cols) +
    theme(plot.title=element_text(size=24, face="bold", hjust=.5),
          axis.title=element_text(size=20),
          axis.text.x=element_text(size=20),
          axis.ticks.x=element_blank(), 
          axis.text.y=element_text(size=20)) +
    labs(title=paste(name, 'chromosome', chr.num))
  print(p)
  dev.off()
}
