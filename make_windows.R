
# Make a bed file of variable sized windows. Each window should contain approximately the same 
# number of reads. 

# Usage: Make_windows.R <reads>.bam <number of reads per window>

# Example: Make_windows.R MAPPED_MASKED/xx-euploid_all.q30.rmdup.bam 1000

library(GenomicAlignments)

args <- commandArgs(trailingOnly = TRUE)
bam.file <- args[1]
n <- as.numeric(args[2])

reads <- readGAlignments(bam.file)
mcols(reads)$mid <- (end(reads)+start(reads)) %/% 2

# Split up reads by chromosome
read.list <- split(reads, seqnames(reads))

# Print the file header
cat(paste0('chr', '\t', 'start', '\t', 'end', '\t', 'count', '\n'))

# For each chromosome, split reads into chunks with n reads (middles) in each
for(i in 1:length(read.list)) {
  # Get total reads for this chromosome
  len <- length(read.list[[i]])
  
  # Sort reads by midpoint
  read.list[[i]] <- read.list[[i]][order(mcols(read.list[[i]])$mid)]

  start <- 1
  bin <- read.list[[i]][1:n]
  end <- max(mcols(bin)$mid)
  in.bin <- which(mcols(read.list[[i]])$mid >= start &
                  mcols(read.list[[i]])$mid <= end)
  bin <- read.list[[i]][in.bin]
  remaining <- read.list[[i]][-in.bin]
  cat(paste0(names(read.list)[[i]], '\t', start, '\t', end, '\t', length(bin), '\n'))

  while(length(remaining) > n) {
    start <- end+1
    bin <- remaining[1:n]
    end <- max(mcols(bin)$mid)
    in.bin <- which(mcols(remaining)$mid >= start &
                      mcols(remaining)$mid <= end)
    bin <- remaining[in.bin]
    remaining <- remaining[-in.bin]
    cat(paste0(names(read.list)[[i]], '\t', start, '\t', end, '\t', length(bin), '\n'))
  }

  start <- end+1
  bin <- remaining
  end <- seqlengths(read.list)[[i]]
  in.bin <- which(mcols(remaining)$mid >= start &
                    mcols(remaining)$mid <= end)
  bin <- remaining[in.bin]
  remaining <- remaining[-in.bin]
  cat(paste0(names(read.list)[[i]], '\t', start, '\t', end, '\t', length(bin), '\n'))
}