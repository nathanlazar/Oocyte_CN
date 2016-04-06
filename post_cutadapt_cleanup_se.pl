#!/usr/bin/perl 

################################################################################
#
#       post_cutadapt_cleanup_se.pl
#
#       This script cleans up a fastq file after it's been
#       trimmed of adapter sequences by cutadapt. 
#
#       Usage:   
#       post_cutadapt_cleanup.pl <fastq_file> <5' trim> <3' trim> > filter_stats.txt
#
#       The script produces a new file named from the input file.
#       The new file is the input file name with '.clean' appended to it.
#
#       If the sequence length of either of the read is less than 25 it is
#       not printed to the output (filtered out).
#
#       If the read fails the Illumina chastity filter it is not printed to 
#       the output. The chastity filter flag is a "Y" if it fails
#       so flags that match this expression ' =~ \:Y\:' are filtered out.
#
#       Larry Wilhelm (modified for single end by Nathan Lazar 12/22/2015)
#                     (nathan dot lazar at gmail dot com)
#       4.12.2012
################################################################################

my $f1 = $ARGV[0];  # fastq read pair 1
my $trim_len_5 = $ARGV[1];
my $trim_len_3 = $ARGV[2];

open(F1,$f1) or die "can't open input file $f1 $!\n";
my $out1 = $f1.".clean";
open(OUT1,">$out1") or die "can't open output file $out1 $!\n";

#my ($n,$len_filtered,$chaste_filtered,$c,$Ns_trimmed);
my $n = 0;
my $len_filtered = 0;
my $chaste_filtered = 0;
my $c = 0;
my $Ns_trimmed = 0;

while(<F1>){
        my $f1_hdr = $_;
        my $seq    = <F1>;
        my $hdr1   = <F1>;
        my $qv1    = <F1>;
        
        chomp $seq;
        chomp $qv1;

        $n++;
        if (length($seq) < 25){
                $len_filtered++;
                next;
        }
        
        if ($seq =~ /N/){
				$Ns_trimmed++;
				next;
		}
        
        

        my ($j,$I1) = split(/\s+/,$f1_hdr);
        if ($I1 =~ /\:Y\:/){
                $chaste_filtered++;
                next;
        }

        $c++;
        
        # Trim 5' ends
        my $s1 = substr($seq,$trim_len_5,length($seq) - $trim_len_5);
        my $q1 = substr($qv1,$trim_len_5,length($qv1)- $trim_len_5);
        
        # Trim 3' ends
        $s1 = substr($s1,0,length($s1) - $trim_len_3);
        $q1 = substr($q1,0,length($q1) - $trim_len_3);

        print OUT1 $f1_hdr.$s1."\n".$hdr1.$q1."\n";
}

print "Total=$n, len_filtered=$len_filtered, Ns_filtered=$Ns_trimmed, chaste_filtered=$chaste_filtered,remaining=$c\n";
