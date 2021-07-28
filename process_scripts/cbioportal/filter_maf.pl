#!/usr/bin/perl -w
#patient_sample_uuid.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'prjid|p=s');

my @maffiles = @ARGV;

open MAFOUT, ">variants.maf" or die $!;
my @mincols = ('Hugo_Symbol','Entrez_Gene_Id','Variant_Classification',
	       'Tumor_Sample_Barcode','HGVSp_Short','Protein_position',
	       'SWISSPROT','t_alt_count','t_ref_count','n_alt_count','n_ref_count');
#($hash{'Hugo_Symbol'},$hash{'Entrez_Gene_Id'},$hash{'Variant_Classification'},$hash{'Tumor_Sample_Barcode'},$hash{'HGVSp_Short'},$hash{'Protein_position'},$hash{'SWISSPROT'});

foreach $maf (@maffiles) {
  open MAF, "<$maf" or die $!;
  while (my $line = <MAF>) {
    chomp($line);
    if ($line =~ m/#/) {
      print MAFOUT $line,"\n";
    }
    elsif ($line =~ m/Hugo_Symbol/i) {
      @mafcols = split(/\t/,$line);
      print MAFOUT join("\t",@mincols),"\n";
    }else {
      my @row = split(/\t/,$line);
      my %hash;
      foreach my $i (0..$#mafcols) {
	$row[$i]  = '' unless $row[$i];
	$hash{$mafcols[$i]} = $row[$i];
      }
      next if ($hash{Variant_Classification} =~ m/Silent|Intron|UTR|Flank|IGR|RNA|Splice_Region/);
      next unless ($hash{FILTER} =~ m/PASS/);
      $mafids{$hash{Tumor_Sample_Barcode}} = 1;
      my @newline;
      foreach $i (0..$#mincols) {
	  push @newline, $hash{$mincols[$i]};
      }
      print MAFOUT join("\t",@newline),"\n";
    }
  }
  close MAF;
}
close MAFOUT;
open SEQD, ">case_lists/sequenced.txt" or die $!;
print SEQD join("\n","cancer_study_identifier: $opt{prjid}",
		"stable_id: $opt{prjid}_sequenced",
		"case_list_name: Sequenced",
		"case_list_description: Sequenced Samples",
		"case_list_ids:".join("\t",keys %mafids)),"\n";
close SEQD;
