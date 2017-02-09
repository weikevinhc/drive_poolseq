# drive_poolseq
Perl and R codes for review.

After the processing with GATK, the merged VCF file which contains all parental samples and crosses will be subject to the perl scripts in the following way.

For expected heterozygous sites from parental lines:
vcf.alt.homo.pl will divide up the the samples into individual vcf files leaving only alternative homozygous sites.
To generate the expected homozygous sites exp.het.pl is used with the vcf files from P1 and P2.

To process the crosses and then polarize with parental the expecte het sites:
vcf.het.parser.pl will divide the samples into individual vcf files leaving only heterozygous sites.
The individual vcf files for the crosses are then processed with cross.het.pl which trims the files into a more maneageable format (.het.txt).

cross.polar.pl will then polarize the .het.txt file with the expected het sites, generating three files .shared.txt, .shared.miss.txt, .shared.ex.txt. The first one contains polarized het sites that are found in both data sets and allele counts at each site; it will then be imported into R. The .shared.miss.txt contains expected het sites missing in the .het.txt and vice versa for .shared.ex.txt.

The R script drive_new.R contains the codes used for trimming the SNPs based on read depth, bias correction, recombination decay estimates, and plotting the data.
