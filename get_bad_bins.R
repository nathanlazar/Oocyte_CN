# nathan dot lazar at gmail dot com

# Determine which bins are consistently outliers and write them
# to std out

# Usage: get_bad_bins.R MAPPED_Batch1/*200_counts.txt > 
#          Mapped_Batch1_bad_bins_200.txt

library(dplyr)    # General data frame manipulation

# Load data
args <- commandArgs(trailingOnly = TRUE)
file.pattern <- args[1]
in.dir <- dirname(file.pattern)
files <- list.files(path = in.dir, 
                    pattern = basename(file.pattern))

# Read in first file to determine size of data frame
f1 <- read.table(file.path(in.dir, files[1]), 
                 header=T, stringsAsFactors=F)

# Make data frame storing bin loci
bins <- dplyr::select(f1, chr, start, end)

# Make matrices to store bin counts and expected percentages
# for all samples
count.mat <- matrix(0, nrow=nrow(f1), ncol=length(files))
exp.per.mat <- count.mat

count.mat[,1] <- f1$count
exp.per.mat[,1] <- f1$expected

# Read in counts and fill count matrix
for(i in 2:length(files)) {
  counts <- read.table(file.path(in.dir, files[i]),
                       header=T, stringsAsFactors=F)
  count.mat[,i] <- counts$count
  exp.per.mat[,i] <- counts$expected
}

# Add one to all counts (to avoid zero effects)
count.mat <- count.mat + 1

# Make a matrix of the percentage of reads in each bin
per.mat <- sweep(count.mat, 2, apply(count.mat, 2, sum), FUN="/")

# Make a matrix of the ratios of observed to expected reads 
# in each bin
ratio.mat <- per.mat/exp.per.mat

# Get the mean ratio for each bin
bins$mean.ratios <- apply(ratio.mat, 1, mean)

# Get bounds (mean of middle 99% of data +/- 3 s.d.) 
# for each chromosome
chroms <- data.frame(chr=unique(bins$chr), 
                     mean=0, sd=0,
                     min=0, max=0)
for(i in 1:nrow(chroms)) {
  all.ratios <- bins$mean.ratios[bins$chr==chroms$chr[i]]
  ratios.99 <- all.ratios[all.ratios > quantile(all.ratios, 0.005) &
                          all.ratios < quantile(all.ratios, 0.995)]
  chroms$mean[i] <- mean(ratios.99)
  chroms$sd[i] <- sd(ratios.99)
  chroms$min[i] <- chroms$mean[i] - 3*chroms$sd[i]
  chroms$max[i] <- chroms$mean[i] + 3*chroms$sd[i]
}

# Mark bins that are more than 2 standard deviations from 
# the mean of the mean ratios for each chromosome as "bad"
bins <- merge(bins, chroms, by='chr', all=T)
bins$bad <- 0
bins$bad[bins$mean.ratios < bins$min] <- 1
bins$bad[bins$mean.ratios > bins$max] <- 1

write.table(dplyr::select(bins, chr, start, end, bad),
            file="", quote=F, row.names=F, sep='\t')