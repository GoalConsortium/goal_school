#!/usr/bin/perl -w
#translocation2cbioportal.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;

my $results= GetOptions (\%opt,'datadir|d=s','gbuilddir|g=s');

open ENT_GENE, "<$opt{datadir}\/gene_info.human.txt" or die $!;
my %entrez;
my %entgene;
my $ent_header = <ENT_GENE>;
while (my $line = <ENT_GENE>){
  chomp $line;
  my @row = split(/\t/, $line);
  $entrez{$row[2]}=$row[1];
}
close ENT_GENE;

my @fusion_files = @ARGV;
open OUT, ">variants.fusion.txt" or die $!;
print OUT join("\t",'Hugo_Symbol','Entrez_Gene_Id','Center',
	       'Tumor_Sample_Barcode','Fusion','DNA_support',
	       'RNA_support','Method','Frame','Fusion_Status');

foreach my $ffile (@fusion_files) {
    open FF, "<$ffile" or die $!;
    my $head = <FF>;
    chomp($head);
    my @colnames = split(/\t/,$head);
    $fname = basename($ffile);
    my $sample = (split(/\./,$fname))[0];
    while (my $line = <FF>) {
	chomp($line);
	my @row = split(/\t/,$line);
	my %hash;
	foreach my $i (0..$#row) {
	    $hash{$colnames[$i]} = $row[$i];
	}
	next unless $hash{Filter} eq 'PASS';
	print OUT join("\t",$hash{LeftGene},$entrez{$hash{LeftGene}},
		       '',$sample,$hash{FusionName},$hash{DNAReads},
		       $hash{RNAReads},'StarFusion',$hash{FusionType},
		       uc($hash{SomaticStatus})),"\n" if $entrez{$hash{LeftGene}};
	print OUT join("\t",$hash{RightGene},$entrez{$hash{RightGene}},
		       '',$sample,$hash{FusionName},$hash{DNAReads},
		       $hash{RNAReads},'StarFusion',$hash{FusionType},
		       uc($hash{SomaticStatus})),"\n" if $entrez{$hash{RightGene}};
	
    }

}
