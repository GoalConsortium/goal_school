#!/usr/bin/perl -w 
#cbioPortal_documents.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;

my $results= GetOptions (\%opt,'fpkm|f=s','logcpm|l=s','cnv|c=s','prefix|p=s','help|h','datadir|d=s');

open ENT_GENE, "<$datadir\/gene_info.human.txt" or die $!;
my %entrez;
my %entgene;
my $ent_header = <ENT_GENE>;
while (my $line = <ENT_GENE>){
  chomp $line;
  my @row = split(/\t/, $line);
  $entgene{$row[6]}{$row[2]}=$row[1];
  #$entrez{$row[2]}=$row[1];
}
close ENT_GENE;
open ENT_ENS, "<$datadir\/genenames.txt" or die $!;
my $gn_header = <ENT_ENS>;
my %ensym;
while (my $line = <ENT_ENS>){
  chomp $line;
  my @row = split(/\t/, $line);
  $entrez{$row[3]}=$entrez{$row[4]};
}
close ENT_ENS;
open ENT_ENS, "<$datadir\/gene2ensembl.human.txt" or die $!;
my $ens_header = <ENT_ENS>;
while (my $line = <ENT_ENS>){
  chomp $line;
  my @row = split(/\t/, $line);
  $entrez{$row[2]}=$row[1];
}
close ENT_ENS;

if($opt{fpkm}){
  open FPKM, "<$opt{fpkm}" or die $!;
  open OUTF, ">$opt{prefix}\.data_fpkm.cbioportal.txt" or die $!;
  print OUTF join("\t","Entrez_Gene_Id",$opt{prefix}),"\n";
  my %fpkm;
  my $fpkm_header = <FPKM>;
  while(my $line = <FPKM>){
    chomp $line;
	my $entrezid=0;
    my ($id,$gene,$chr,$strand,$start,$end,$coverage,$fpkm,$tpm) = split(/\t/,$line);
	$chr =~ s/^chr//g;
    my $ensembl = (split(/\./,$id))[0];
    if ($entrez{$ensembl}) {
      $entrezid = $entrez{$ensembl};
    }elsif($entgene{$chr}{$gene}){
      $entrezid = $entgene{$chr}{$gene};
    }
    next unless ($entrezid);
    print OUTF join("\t",$entrezid,$fpkm),"\n"; 
  }
  close OUTF;
}

if($opt{logcpm}){
  open IN, "<$opt{logcpm}" or die $!;
  open OUTL, ">$opt{prefix}\.data_logCPM.cbioportal.txt" or die $!;
  print OUTL join("\t","Entrez_Gene_Id",$opt{prefix}),"\n";
  $fname = basename($opt{logcpm});
  my $sample = (split(/\./,$fname))[0];
  my $command = <IN>;
  my $head = <IN>;
  chomp($head);
  my $total = 0;
  while (my $line = <IN>) {
    chomp($line);
    my @row = split(/\t/,$line);
    my $gene = $row[0];
	my $chrom = (split(";",$row[1]))[0];
	$chrom =~ s/^chr//g;
    my $ct = $row[-1];
    next if($gene =~ m/^__/);
    $cts{$chrom}{$gene} = $ct;
    $total += $ct;
  }
	print $total."\n";
  close IN;
  foreach my $ens_chr (keys %cts) {
	foreach my $ens (keys $cts{$ens_chr}){
    	next unless $entgene{$ens_chr}{$ens};#($entrez{$ens} or $entgene{$gene}{$chrom});
    	#unless ($cts{$ens_chr}{$ens}){
    	#  $cts{$ens_chr}{$ens} = 0;
    	#}
    	$cpm = ($cts{$ens_chr}{$ens}/$total)*1e6;
    	print OUTL join("\t",$entgene{$ens_chr}{$ens},sprintf("%.2f",log2($cpm))),"\n";
	}
  }
  close OUTL;
}

sub log2 {
    $n = shift @_;
    if ($n < 1) {
	return 0;
    }else {
	return(log($n)/log(2));
    }
}
