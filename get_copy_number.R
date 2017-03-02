# nathan dot lazar at gmail dot com

# Uses CBS algorithm and package DNAcopy to estimate copy numbers given 
# ratios calculated with count_bin_reads.R. Copy number estimates are 
# an added column in the given matrix of counts and ratios.

# Usage: Rscript get_copy_number.R <bin_count_file.counts> 

# Example: Rscript get_copy_number.R \
#          COUNTS/500/rh150409-1-b1-B1_S1.counts

library(dplyr)    # General data frame manipulation
library(DNAcopy)  # Implements circular binary segmentation


# Functions ---------------------------------------------------------------

get_seg <- function(counts, base.cn) {
  # Use circular binary segmentation to call copy number states 
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

  counts.seg$output$round <- round(counts.seg$output$seg.mean)
  counts.seg
}

# Recalibrate normalization based on the segmentation counts --------------

recalibrate <- function(counts, counts.seg, bin.size=200, base.cn=2, verbose=T) {
  # Assign putative copy number to each window
  counts$cn <- unlist(mapply(rep, counts.seg$round, counts.seg$num.mark))
  
  if(verbose) {
    print('Putative copy number means before correction:')
    counts %>% tbl_df %>% group_by(cn) %>% 
      summarise(mean.ratio=mean(ratio)) %>% print
  }
  
  # Correct for the differing amount of DNA in aneuploid samples
  control.tot <- sum(counts$count!=0) * bin.size * base.cn
  samp.tot <- sum(counts$cn[counts$count!=0])* bin.size/base.cn * 2
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

if(grepl('500', count.file)) bin.size=500
if(grepl('1000', count.file)) bin.size=1000
if(grepl('2000', count.file)) bin.size=2000
if(grepl('4000', count.file)) bin.size=4000

counts <- read.table(count.file, header=T, stringsAsFactors=F)

# Remove small chromosomes (< 50 bins)
chr.bins <- counts %>% tbl_df %>%
  group_by(chr) %>%
  summarise(bins=length(count))

counts <- counts[counts$chr %in% chr.bins$chr[chr.bins$bins >=50],]
chrs <- unique(counts$chr)
chr.sort <- c(paste0('chr', sort(as.numeric(sub('chr', '', chrs[chrs!='chrX'])))),
            'chrX')

counts$chr <- factor(counts$chr, levels=chr.sort)

out <- sub('.counts', '.cn', count.file)
name <- sub('.counts', '', count.file)
split.name <- strsplit(name, split='/', fixed=T)[[1]]
name <- split.name[length(split.name)]

# mean(counts$ratio)

# Remove M and Y chroms
tot.reads <- sum(counts$count)
reads.rem <- sum(counts$count[counts$chr %in% c('chrM', 'chrY')])
per.rem <- reads.rem/tot.reads
exp.per.rem <- sum(counts$expected[counts$chr %in% c('chrM', 'chrY')])
print(sprintf('%.2f%% of reads removed that mapped to chromosomes M and Y', 
              per.rem * 100))
adj.fact <- (1-exp.per.rem)/(1-per.rem)
print(sprintf('Expected %.2f%% of reads to map to chromosomes M and Y so we correct by a factor of %.2f', 
              exp.per.rem * 100, adj.fact))
counts <- counts[!(counts$chr %in% c('chrM', 'chrY')),]
counts <- droplevels(counts)

# Adjust and ratios for removed reads
counts$per <- counts$per * adj.fact
counts$ratio <- counts$per/counts$expected

# Reorder chromosomes
counts <- with(counts, counts[order(chr),])

# Set bins w/ less than 2 reads to zero
read.cut <- 2
per.rem <- sum(counts$count[counts$count < read.cut])/sum(counts$count)
print(sprintf('%.2f%% of reads removed in bins with less than %s reads', 
              per.rem * 100, read.cut))
counts$count[counts$count < read.cut] <- 0
# Adjust expected counts removing bins w/ zero counts
counts$expected <- counts$expected * (1-per.rem/sum(counts$count))
counts$ratio[counts$count < read.cut] <- 0
counts$per[counts$count < read.cut] <- 0

# If most bins are empty we expect a base copy number of 1 
# otherwise base copy number is set to 2
if(mean(counts$count == 0) > .5) {
  base.cn <- 1
} else base.cn <- 2

# Get segmentation for each chromosome separately
chroms <- levels(counts$chr)
seg.list <- list()
for(chrom in chroms) {
  count.chr <- filter(counts, chr==chrom)
  seg.list[[chrom]] <- get_seg(count.chr, base.cn)$output
}
counts.seg <- seg.list[[1]]
for(i in 2:length(seg.list))
  counts.seg <- rbind(counts.seg, seg.list[[i]])

# Recalibrate counts given the putative copy numbers
counts <- recalibrate(counts, counts.seg, bin.size, base.cn, verbose=T)

# Re-call copy numbers?
# seg.list2 <- list()
# for(i in 1:nrow(bins.per.chrom)) {
#   chrom <- bins.per.chrom$chr[i]
#   count.chr <- filter(counts, chr==chrom)
#   seg.list2[[chrom]] <- get_seg(count.chr, 1)$output
# }
# counts.seg2 <- seg.list2[[1]]
# for(i in 2:length(seg.list2))
#   counts.seg2 <- rbind(counts.seg2, seg.list2[[i]])
# 
# plot(counts.seg2$round)

# Renumber the segmentation x-values
for(i in 2:nrow(counts.seg)) {
  counts.seg$loc.start[i] <- counts.seg$loc.end[i-1]+1
  counts.seg$loc.end[i] <- counts.seg$loc.start[i] + counts.seg$num.mark[i] - 1
}

# Add index for whole genome and each chrom
counts$idx <- 1:nrow(counts)
bins.per.chrom <- counts %>% tbl_df %>% group_by(chr) %>% summarise(len=length(chr)) 
counts$bin <- unlist(sapply(bins.per.chrom$len, function(x) 1:x))

# Add copy number into counts
counts$cn <- 0
for(i in 1:nrow(counts.seg)) 
  counts$cn[counts$idx >= counts.seg$loc.start[i] & 
            counts$idx <= counts.seg$loc.end[i]] <- counts.seg$round[i]

# Save count data.frame
write.table(counts, file=out)
