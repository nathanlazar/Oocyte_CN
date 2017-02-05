#!/bin/bash

bin=/opt/installed
bin2=/home/exacloud/lustre1/users/lazar/bin/bedtools2/bin
cores=12
REF=/home/groups/hoolock/u1/nl/Oocyte_CN/RheMac8/rheMac8.fa
out=MAPPED_SE
F1=$1

# This script runs bwa-mem and filters for quality of 30 (which excludes multi-mapping reads)
# Results are placed in MAPPED_SE
# Example outputs are:
#  reads.bam
#  reads.q30.sort.bam
#  reads.q30.sort.bam.bai
#  reads.mapped_reads   # Read counts at each step and summary measures

NM=${F1/.trim.uniq.fq.gz/}
NM=${NM/_L001/}
NM=$(basename $NM)

# Get basic stats from the fastq file
echo 'total unique per_unique most_common_seq most_com_seq_count per_most_common' > $out/$NM.mapped_reads
zcat $F1 | awk '(NR%4==0){read=$1;total++;count[read]++} \
  END{for(read in count){if(!max||count[read]>max) \
    {max=count[read];maxRead=read};\
    if(count[read]==1){unique++}}; \
  print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}' >> $out/$NM.mapped_reads

# Map reads, filter out secondary alignments and unmapped reads
$bin/bwa mem -M -t $cores $REF $F1 > $out/$NM.sam 
$bin/samtools view -bShF0x0104 $out/$NM.sam > $out/$NM.bam
# -M marks shorter split hits as secondary
# -t is the number of threads
# -b output bam
# -S input is sam
# -h print sam header
# -F 0x0100 skip secondary hits 0x0004 skip unmapped 

# Record number of mapped reads
echo -n "Uniquely mapped reads: " >> $out/$NM.mapped_reads
$bin/samtools view $out/$NM.bam | wc -l >>  $out/$NM.mapped_reads

# Filter by mapping quality scores over 30
$bin/samtools view -b -q 30 $out/$NM.bam  > $out/$NM.q30.bam
echo -n "Passing Q 30 filter: " >> $out/$NM.mapped_reads
$bin/samtools view $out/$NM.q30.bam | wc -l >> $out/$NM.mapped_reads

# Sort
$bin/samtools sort -@ $cores $out/$NM.q30.bam $out/$NM.q30.sort

# Index the filtered and sorted .bam file
$bin/samtools index $out/$NM.q30.sort.bam

echo "#########################" >> $out/$NM.mapped_reads

# Count mapping coordinates, unique mapping coordinates & maximum number of reads
# mapped to the same position
echo "Mapping positions:" >> $out/$NM.mapped_reads
echo "Total Unique PerUnique MaxCoor CountMaxCoor PerMaxCoor" >> $out/$NM.mapped_reads
$bin2/bamToBed -i $out/$NM.q30.sort.bam | \
awk '{coordinates=$1":"$2"-"$3; total++; count[coordinates]++}\
  END{for(coordinates in count){if(!max||count[coordinates]>max) \
    {max=count[coordinates];maxCoor=coordinates}; \
    if(count[coordinates]==1){unique++}}; \
  print total,unique,unique*100/total,maxCoor,count[maxCoor],count[maxCoor]*100/total}' >> $out/$NM.mapped_reads

# Clean up intermediate files
rm $out/$NM.sam
rm $out/$NM.bam
rm $out/$NM.q30.bam

