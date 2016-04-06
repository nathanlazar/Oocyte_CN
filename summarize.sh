#!/bin/bash

cd MAPPED

grep -A1 total *.mapped_reads | \
  awk '((NR-2)%3==0) {print $1}' | \
  sed 's/.mapped_reads-/ /' | sort | uniq > read_stats.txt

echo "#################################" >> read_stats.txt
echo "Uniquely mapped reads" >> read_stats.txt

grep  Uniquely *.mapped_reads | \
  sed 's/.mapped_reads:Uniquely mapped reads://' | sort | uniq >> read_stats.txt

echo "#################################" >> read_stats.txt
echo "Passing Q30 filter" >> read_stats.txt

grep  Passing *.mapped_reads | \
  sed 's/.mapped_reads:Passing Q 30 filter://' | sort | uniq >> read_stats.txt

echo "#################################" >> read_stats.txt
echo "Unique mapping positions" >> read_stats.txt

grep -A1 Total *.mapped_reads | \
  awk '((NR-2)%3==0) {print $1,$2}' | \
  sed -r 's/.mapped_reads-[0-9]+/ /' | sort | uniq >> read_stats.txt
