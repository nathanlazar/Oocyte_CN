#!/bin/bash

REF=$1
#example: /u0/dbase/genomes/rhesus_v7/MaSuRCA_Rhesus_Genome_v7_20130927.fasta
F1=$2
F2=$3
NM=$4

# Get basic stats from the fastq file
echo 'total unique per_unique most_common_seq most_com_seq_count per_most_common' > $NM.mapped_reads
zcat $F1 | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}' >> $NM.mapped_reads
zcat $F2 | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}' >> $NM.mapped_reads

# Map reads
bwa mem -t 8 $REF $F1 $F2 | samtools view -bSh -  > $NM".bam"

# Sort
samtools sort $NM".bam" $NM

# Index
samtools index $NM".bam"
echo "All reads:" > $NM".mapped_reads"
samtools view $NM".bam" | wc -l >>  $NM".mapped_reads"

# Filter by quality scores over 30
samtools view -b -q 30 $NM".bam"  > $NM".q30.bam" 
echo "Passing Q 30 filter:" >> $NM".mapped_reads"
samtools view $NM".q30.bam" | wc -l >> $NM".mapped_reads"

# Remove duplicates w/ samtools rmdup -S
# This treats reads as single-ended which eliminates big read pileups
/usr/local/bin/samtools-0.1.9/samtools rmdup -S $NM".q30.bam"  $NM".q30.rmdup.bam"
echo "After samtools rmdup:" >> $NM".mapped_reads"
samtools view $NM".q30.rmdup.bam" | wc -l >> $NM".mapped_reads"

