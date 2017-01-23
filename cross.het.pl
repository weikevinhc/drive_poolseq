#!usr/bin/perl
# cross.het.pl
# determine and filter for heterozygous SNPs from cross.vcf. 
use strict;
use warnings;
use POSIX;

my $input1 = shift @ARGV; # vcf file
my $input2 = shift @ARGV; # output base name
my $GQ = 20; ## phred quality score

open VCF1, $input1 or die;
open OUTPUT, ">$input2.het.txt" or die;
print OUTPUT "chr\tpos\tref\talt\tGT\n";	
while (<VCF1>) {
	next if $_ =~ /^#/;
	s/[\r\n]+$//;
	my @line = split("\t", $_);
	next if (length($line[4]) > 1);
	next if (length($line[3]) > 1);
	my @format = split(":", $line[8]);
	unless ($format[3]) {
		next;
	}
	next if ($format[3] ne "GQ");
	my @GT = split(":", $line[9]);
	my @counts = split(",", $GT[1]);
	next if (scalar(@counts) < 2);
	next if ($counts[0] == 0 || $counts[1] == 0);
	if ($GT[0] eq "0\/1") {
		if ($GT[3] > $GQ) {
			print OUTPUT "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$line[9]\n";
		}		
	}	
}

	
			
