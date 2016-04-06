#!/usr/bin/perl 

################################################################################
#
#       post_cutadapt_cleanup.pl
#
#       This script cleans up a pair of fastq files after they've been
#       trimmed of adapter sequences by cutadapt. The script is designed to only 
#       work on properly synchronized pairs of input fastq files. 
#
#       Usage:   
#       post_cutadapt_cleanup.pl <fastq_pair_1> <fastq_pair_2> > filter_stats.txt
#
#       The script produces two new files that are named from the input files.
#       The new files are the input file name with '.clean' appended to it.
#
#       If the sequence length of either of the read pairs is less than 25 both
#       members of the pair are not printed to the output (filtered out).
#
#       If either read pair fails the Illumina chastity filter both members are
#       not printed to the output. The chastity filter flag is a "Y" if it fails
#       so flags that match this expression ' =~ \:Y\:' are filtered out.
#
#       The output should therefore be two properly synchronized fastq files of 
#       read pairs.
#
#       Larry Wilhelm
#       4.12.2012
################################################################################

my $f1 = $ARGV[0];  # fastq read pair 1
my $f2 = $ARGV[1];  # fastq read pair 2
my $trim_len_5 = $ARGV[2];
my $trim_len_3 = $ARGV[3];
my $trim_len_5_2 = $ARGV[4];
my $trim_len_3_2 = $ARGV[5];


open(F1,$f1) or die "can't open input file $f1 $!\n";
open(F2,$f2) or die "can't open input file $f2 $!\n";
my $out1 = $f1.".clean";
my $out2 = $f2.".clean";
open(OUT1,">$out1") or die "can't open output file $out1 $!\n";
open(OUT2,">$out2") or die "can't open output file $out2 $!\n";

my ($n,$len_filtered,$chaste_filtered,$c,$Ns_trimmed);

while(<F1>){
        my $f1_hdr = $_;
        my $seq    = <F1>;
        my $hdr1   = <F1>;
        my $qv1    = <F1>;

        my $f2_hdr = <F2>;
        my $seq2   = <F2>;
        my $hdr2   = <F2>;
        my $qv2    = <F2>;
        
        chomp $seq;
        chomp $qv1;
        chomp $seq2;
        chomp $qv2;

        $n++;
        if (length($seq) < 25 || length($seq2) < 25){
                $len_filtered++;
                next;
        }
        
        if ($seq =~ /N/ || $seq2 =~ /N/){
				$Ns_trimmed++;
				next;
		}
        
        

        my ($j,$I1) = split(/\s+/,$f1_hdr);
        my ($k,$I2) = split(/\s+/,$f2_hdr);
        if ($I1 =~ /\:Y\:/ || $I2 =~ /\:Y\:/){
                $chaste_filtered++;
                next;
        }

        $c++;
        
        # Trim 5' ends
        my $s1 = substr($seq,$trim_len_5,length($seq) - $trim_len_5);
        my $q1 = substr($qv1,$trim_len_5,length($qv1)- $trim_len_5);
        my $s2 = substr($seq2,$trim_len_5_2,length($seq2) - $trim_len_5_2);
        my $q2 = substr($qv2,$trim_len_5_2,length($qv2)- $trim_len_5_2);
        
        # Trim 3' ends
        $s1 = substr($s1,0,length($s1) - $trim_len_3);
		$q1 = substr($q1,0,length($q1) - $trim_len_3);
		$s2 = substr($s2,0,length($s2) - $trim_len_3_2);
		$q2 = substr($q2,0,length($q2) - $trim_len_3_2);

        print OUT1 $f1_hdr.$s1."\n".$hdr1.$q1."\n";
        print OUT2 $f2_hdr.$s2."\n".$hdr2.$q2."\n"; 
}

print "Total=$n, len_filtered=$len_filtered, Ns_filtered=$Ns_trimmed, chaste_filtered=$chaste_filtered,remaining=$c\n";
