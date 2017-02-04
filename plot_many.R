# nathan dot lazar at gmail dot com

# Plots the ratio of observed over expected read counts for the whole genome
# for all the files in a directory and combine them with 16 on a page

# Usage: Rscript plot_many.R <dir>

# Example: Rscript plot_many.R ''D:/Box Sync/Rhesus_Embryos/Oocyte_CN/PE/1000/'

library(dplyr)    # data frame manipulation
library(ggplot2)  # Used for plotting
library(grid)     # Used for plotting 

# Plot function ----------------------------------------------------------------

plot_wg <- function(cn.file) {
  name <- sub('.cn', '', cn.file)
  split.name <- strsplit(name, split='/', fixed=T)[[1]]
  name <- split.name[length(split.name)]

  if(grepl('500', cn.file)) bin.size=500
  if(grepl('1000', cn.file)) bin.size=1000
  if(grepl('2000', cn.file)) bin.size=2000
  if(grepl('4000', cn.file)) bin.size=4000
  
  counts <- read.table(cn.file, header=T, stringsAsFactors=F)
  
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
  counts$chr <- factor(counts$chr, levels=c(paste0('chr', 1:20), 'chrX'))
  counts.seg$chr <- factor(counts.seg$chr, levels=levels(counts$chr))
  counts.seg$cn <- as.numeric(counts.seg$cn)
  counts.seg$start <- as.numeric(counts.seg$start)
  counts.seg$end <- as.numeric(counts.seg$end)
  
  # Prep for plotting
  bins.per.chrom <- counts %>% tbl_df %>% group_by(chr) %>% summarise(len=length(chr)) 
  bins.per.chrom$mid <- cumsum(bins.per.chrom$len) - (bins.per.chrom$len/2)
  n.chr <- length(unique(counts$chr))
  counts$col <- 2
  counts$col[counts$chr %in% unique(counts$chr)[seq(1,n.chr,2)]] <- 1
  counts$col <- as.factor(counts$col)
  ylim <- c(0, min(max(counts$ratio), 6))
  
  # Plot whole genome
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

  return(p)
  # gt <- ggplot_gtable(ggplot_build(p))
  # gt$layout$clip[gt$layout$name == "panel"] <- "off"
  # grid.draw(gt)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Main script -------------------------------------------------------------

# Load data
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]

files <- list.files(dir, '*.cn')

k <- 1
for(i in 1:(length(files) %/% 16)) {
  png(file=paste0(dir, '/multi_', i, '.png'), width=480*8, height=480*4)
  plots <- list()
  for(j in 1:16) {
    file.cn <- paste0(dir, files[k])
    plots[[j]] <- plot_wg(file.cn)
    k <- k+1
  }
  multiplot(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
            plots[[5]], plots[[6]], plots[[7]], plots[[8]], 
            plots[[9]], plots[[10]], plots[[11]], plots[[12]], 
            plots[[13]], plots[[14]], plots[[15]], plots[[16]], cols=4)
  dev.off()
}
