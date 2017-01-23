#!usr/bin/perl
# vcf.het.parser.pl
# separate merged VCF files into individual files for each cross while only retaining het sites. The output names will be based on the sample names in VCF with .het.vcf.
use strict;
use warnings;
use POSIX;


my $input1 = shift @ARGV;
open VCF, $input1 or die;
my @header = ();
while (<VCF>) {
	next if $_ =~ /^##/;
	s/[\r\n]+$//;
	my @line = split("\t", $_);
	if ($_ =~ /^#/) {
		foreach my $i (9..(scalar(@line)-1)) {
			@header = @line;
			open OUTPUT, ">$header[$i].het.vcf";
			print OUTPUT $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", $line[4], "\t", $line[5], "\t", $line[6], "\t", $line[7], "\t", $line[8], "\t", "GENOTYPE", "\n";
			close OUTPUT;
		}
		next;
	}
	foreach my $i (9..(scalar(@line)-1)) {
		my @GT = split(":", $line[$i]);
		if ($GT[0] eq "0\/1") {
			open OUTPUT, ">>$header[$i].het.vcf";
			print OUTPUT $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", $line[4], "\t", $line[5], "\t", $line[6], "\t", $line[7], "\t", $line[8], "\t", $line[$i], "\n";
		}
	}
}
			
			
			
