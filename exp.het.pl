#!usr/bin/perl
# exp.het.pl
# determine the expected heterozygous sites based on G0 cross. 
use strict;
use warnings;
use POSIX;

my $input1 = shift @ARGV; 
my $input2 = shift @ARGV;
my $input3 = shift @ARGV;
my $GQ = 23; ## phred quality score for genotype call
print "Enter P1 name:\t";
my $strainA = <STDIN>;
chomp $strainA;
print "Enter P2 name:\t";
my $strainB = <STDIN>;
chomp $strainB;

open VCF1, $input1 or die;

my @header = ();
my %pos1 = ();
my %chr_hash = ();
my @chr_array;
warn "opening vcf1\n";
while (<VCF1>) {
	next if $_ =~ /^#/;
	s/[\r\n]+$//;
	my @line = split("\t", $_);
	unless ($chr_hash{$line[0]}) {
		$chr_hash{$line[0]} = 1;
		push @chr_array, $line[0];
	}
	next if (length($line[4]) > 1);
	next if (length($line[3]) > 1);
	my @format = split(":", $line[8]);
	unless ($format[3]) {
		next;
	}
	next if ($format[3] ne "GQ");
	my @GT = split(":", $line[9]);
	if ($GT[0] eq "1\/1") {
		my $value = join(";", $line[3], "1", $line[4], $line[9]);
		$pos1{$line[0]}{$line[1]} = $value;
	}	
}
close VCF1;

open VCF2, $input2 or die;
warn "opening vcf2\n";
while (<VCF2>) {
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
	if ($GT[0] eq "1\/1") {
		my $value = join(";", $line[3], "2", $line[4], $line[9]);
		if ($pos1{$line[0]}{$line[1]}) {
			my @array = split(";", $pos1{$line[0]}{$line[1]});
			if ($array[2] eq $line[4]) {
				delete $pos1{$line[0]}{$line[1]};
			} else {
				$pos1{$line[0]}{$line[1]} = join(";", $pos1{$line[0]}{$line[1]}, $value);
			}
		} else {
			$pos1{$line[0]}{$line[1]} = $value;
		}		
	}	
}

open OUTPUT, ">$input3" or die;
print OUTPUT "chr\tpos\tref\t$strainA\t$strainB\tGT_$strainA\tGT_$strainB\n";

foreach my $chr (@chr_array) {
	foreach my $snp (sort {$a <=> $b} keys %{$pos1{$chr}}) {
		my @value = split(";", $pos1{$chr}{$snp});
		my @GT = split(":", $value[3]);
		next if ($GT[3] < $GQ);
		if ($value[4]) {
			print "hit";
			print OUTPUT "$chr\t$snp\t$value[0]\t$value[2]\t$value[6]\t$value[3]\t$value[7]\n";
		} else {
			if ($value[1] == 1) {
				print OUTPUT "$chr\t$snp\t$value[0]\t$value[2]\t$value[0]\t$value[3]\t.\n";
			} elsif ($value[1] == 2) {
				print OUTPUT "$chr\t$snp\t$value[0]\t$value[0]\t$value[2]\t.\t$value[3]\n";
			}
		}
	}
}
				
			