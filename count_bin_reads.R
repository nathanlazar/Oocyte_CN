# nathan dot lazar at gmail dot com

# Counts the number of reads in each of supplied bins given a bam file.
# A read is considered in a bin if it's midpoint (or 5' of mid point if even length)
# is in the bin.
# Output is the bin loci, raw counts per bin, percentage of reads in each bin,
# expected percentage of reads in each bin (according to the bin file) and the ratio 
# of the observed percentage of reads in each bin to the expected percentage of 
# reads in each bin. This last value is what is plotted to determine loss of 
# chromosomes (or pieces thereof).

# Usage: Rscript count_bin_reads.R <bin_file.txt> <reads.bam> <output.txt>

# Example: Rscript count_bin_reads.R \
#            xx-bins_200.txt \
#            MAPPED_MASKED/10-24-14-B1_S7_L001.q30.rmdup.bam \
#            MAPPED_MASKED/10-24-14-B1_S7_L001.bin_200_counts.txt

library(GenomicAlignments)

args <- commandArgs(trailingOnly = TRUE)
bin.file <- args[1]
bam.file <- args[2]
out.file <- args[3]

reads <- readGAlignments(bam.file)
# Add midpoints of reads to read object
mcols(reads)$mid <- (end(reads)+start(reads)) %/% 2

bins <- read.table(bin.file, stringsAsFactors=F, header=T, sep='\t')

# Add the percentage of total reads for each bin
bins$per <- bins$count/sum(bins$count)

# Get counts for the new sample
counts <- bins
counts$count <- apply(bins, 1, function(x) length(reads[seqnames(reads) == x[1] &
                  mcols(reads)$mid >= as.numeric(x[2]) &
                  mcols(reads)$mid <= as.numeric(x[3])]))

counts$per <- counts$count/sum(counts$count)
counts$expected <- bins$per
counts$ratio <- counts$per/bins$per

write.table(counts, file=out.file, col.names=T, row.names=F, sep='\t', quote=F)