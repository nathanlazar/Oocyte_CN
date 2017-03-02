# nathan dot lazar at gmail dot com

# Plots the ratio of observed over expected read counts in each of the 
# supplied bins for the whole sample as well as individual chromosomes
# are generated and placed in a directory with the base name of the 
# file plus _plots

# Usage: Rscript plot_copy_number.R <bin_count_file.cn> <out_dir>

# Example: Rscript plot_copy_number.R \
#          COUNTS/500/Monosomy-0_S16.cn'

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

if(grepl('500', cn.file)) bin.size=500
if(grepl('1000', cn.file)) bin.size=1000
if(grepl('2000', cn.file)) bin.size=2000
if(grepl('4000', cn.file)) bin.size=4000

counts <- read.table(cn.file, header=T, stringsAsFactors=F)
                     
# Make output directory if necessary 
if (!file.exists(out.dir))
  dir.create(file.path(out.dir))

# Make a dataframe of the copy number calls
# counts.seg <- tbl_df(counts) %>% 
#   select(chr, idx, cn) %>%
#   group_by(chr, cn) %>%
#   summarise(start=min(idx), end=max(idx))

# Loop to make the counts.seg object. There's probably a way better way to do this
# with dplyr
counts.seg <- data.frame(chr=counts$chr[1], cn=counts$cn[1], 
                         start=counts$idx[1], end=counts$idx[1],
                         stringsAsFactors = F)
j <- 1
for(i in 2:nrow(counts)) {
  chr <- counts$chr[i]
  start <- counts$idx[i]
  end <- counts$idx[i]
  cn <- counts$cn[i]
  if(chr != counts.seg$chr[j]) {  # If new chrom
    j <- j+1
    counts.seg <- rbind(counts.seg, c(chr, cn, start, end))
  } else {
    if(cn == counts.seg$cn[j]) { # If not a new bin
    counts.seg$end[j] <- counts$idx[i]
    } else { # If a new bin
      j <- j+1
      counts.seg <- rbind(counts.seg, c(chr, cn, start, end))
    }
  }
}

chrs <- unique(counts$chr)
chr.sort <- c(paste0('chr', sort(as.numeric(sub('chr', '', 
  chrs[chrs!='chrX'])))), 'chrX')
counts$chr <- factor(counts$chr, levels=chr.sort)
counts.seg$chr <- factor(counts.seg$chr, levels=chr.sort)
counts.seg$cn <- as.numeric(counts.seg$cn)
counts.seg$start <- as.numeric(counts.seg$start)
counts.seg$end <- as.numeric(counts.seg$end)

# Prep for plotting
bins.per.chrom <- counts %>% tbl_df %>% group_by(chr) %>% summarise(len=length(chr)) 
bins.per.chrom$mid <- cumsum(bins.per.chrom$len) - (bins.per.chrom$len/2)
n.chr <- length(chrs)
counts$col <- 2
counts$col[counts$chr %in% chrs[seq(1,n.chr,2)]] <- 1
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
  theme(plot.title=element_text(size=24, face="bold", hjust=.5),
        axis.title=element_text(size=20),
        axis.text.x=element_text(size=14),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_text(size=20)) +
  geom_text(data=bins.per.chrom, size=6, colour='#7F7F7F',
            aes(x=bins.per.chrom$mid, vjust=2.5,
                y=rep(0,nrow(bins.per.chrom)), 
                label=sub('chr', '', bins.per.chrom$chr))) +
  labs(title=name) +
  geom_segment(data=counts.seg, aes(x=start, xend=end, y=cn, yend=cn), size=2)

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)
dev.off()

for(chr in chrs) {
  chr.counts <- counts[counts$chr==chr,]
  chr.counts$idx <- chr.counts$bin
  chr.counts.seg <- counts.seg
  chr.counts.seg <- counts.seg[counts.seg$chr==chr,]
  chr.start <- min(chr.counts.seg$start)
  chr.counts.seg$start <- chr.counts.seg$start - chr.start
  chr.counts.seg$end <- chr.counts.seg$end - chr.start
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
    theme(plot.title=element_text(size=24, face="bold", hjust=.5),
          axis.title=element_text(size=20),
          axis.text.x=element_text(size=20),
          axis.ticks.x=element_blank(), 
          axis.text.y=element_text(size=20)) +
    labs(title=paste(name, 'chromosome', chr.num)) +
    geom_segment(data=chr.counts.seg, aes(x=start, xend=end, y=cn, yend=cn), size=2)
  print(p)
  dev.off()
}
