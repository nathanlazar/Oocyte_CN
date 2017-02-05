
# Make a bed file of variable sized bins. Each bin should contain approximately the same 
# number of reads. 

# Usage: make_bins.R <reads>.bam <number of reads per bin>

# Example: make_bins.R MAPPED_MASKED/xx-euploid_all.q30.rmdup.bam 1000

library(GenomicAlignments)

args <- commandArgs(trailingOnly = TRUE)
bam.file <- args[1]
n <- as.numeric(args[2])

reads <- readGAlignments(bam.file)
mcols(reads)$mid <- (end(reads)+start(reads)) %/% 2

# Subset to just the major chromosomes (without 'Un')
big.chr <- levels(seqnames(reads))
big.chr <- big.chr[!grepl('Un', big.chr)]
reads <- reads[seqnames(reads) %in% big.chr]

# Split up reads by chromosome
read.list <- split(reads, as.character(seqnames(reads)))

# Print the file header
cat(paste0('chr', '\t', 'start', '\t', 'end', '\t', 'count', '\n'))

# For each chromosome, split reads into chunks with n reads (middles) in each
for(i in 1:length(read.list)) {
  
  # Sort reads by midpoint
  read.list[[i]] <- read.list[[i]][order(mcols(read.list[[i]])$mid)]

  end <- 0
  remaining <- read.list[[i]]

  while(length(remaining) > n) {
    start <- end+1
    bin <- remaining[1:n]
    end <- max(mcols(bin)$mid)
    in.bin <- which(mcols(remaining)$mid >= start &
                      mcols(remaining)$mid <= end)
    bin <- remaining[in.bin]
    cat(paste0(names(read.list)[[i]], '\t', start, '\t', end, '\t', length(bin), '\n'))
    remaining <- remaining[-in.bin]
  }

  start <- end+1
  end <- seqlengths(read.list)[[i]]
  bin <- remaining
  cat(paste0(names(read.list)[[i]], '\t', start, '\t', end, '\t', length(bin), '\n'))
}