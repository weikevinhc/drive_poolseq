#!usr/bin/perl
# vcf.snp.parser.pl
# separate merged VCF files into individual files for each parental line while only retaining homozygous alternative sites.
use strict;
use warnings;
use POSIX;


my $input1 = shift @ARGV;
my $GQ = shift @ARGV; genotype quality phred score cutoff
$GQ = 20 unless ($GQ);
open VCF, $input1 or die;
my @header = ();
while (<VCF>) {
	next if $_ =~ /^##/;
	s/[\r\n]+$//;
	my @line = split("\t", $_);
	if ($_ =~ /^#/) {
		foreach my $i (9..(scalar(@line)-1)) {
			@header = @line;
			open OUTPUT, ">$header[$i].alt.vcf";
			print OUTPUT $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", $line[4], "\t", $line[5], "\t", $line[6], "\t", $line[7], "\t", $line[8], "\t", "GENOTYPE", "\n";
			close OUTPUT;
		}
		next;
	}
	foreach my $i (9..(scalar(@line)-1)) {
		my @GT = split(":", $line[$i]);
		if ($GT[0] eq "1\/1" && $GT[3] >= $GQ) {
			open OUTPUT, ">>$header[$i].alt.vcf";
			print OUTPUT $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", $line[4], "\t", $line[5], "\t", $line[6], "\t", $line[7], "\t", $line[8], "\t", $line[$i], "\n";
		}
	}
}
