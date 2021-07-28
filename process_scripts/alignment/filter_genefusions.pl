#!/usr/bin/perl -w
#svvcf2bed.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'fusion|f=s','prefix|p=s','help|h','datadir|r=s');

open OM, "<$opt{datadir}/known_genefusions.txt" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $known{$line} = 1;
}
close OM;

open OM, "<$opt{datadir}/panelgenes.txt" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $keep{$line} = 1;
}

my @exonfiles = `ls */*.exons.txt`;
foreach $efile (@exonfiles) {
    chomp($efile);
    my @leftexons;
    my @rightexons;
    my ($dir,$pair,@etc) = split(/\/|\./,$efile);
    open EFILE, "<$efile" or die $!;
    my $header = <EFILE>;
    while (my $line = <EFILE>) {
	my ($leftgene,$rightgene,$lefttrx,$righttrx,$exonsrc,
	    $exonnum,$exon_chr,$exon_start,$exon_end) = split(/\t/,$line);
	if ($exonsrc =~ m/5/) {
	    push @leftexons, $exonnum;
	}else {
	    push @rightexons, $exonnum;
	}
    }
    if ($leftexons[0] eq  $leftexons[-1]) {
	$exonnuminfo{$dir}{leftgene} = $leftexons[0];
    }else {
	$exonnuminfo{$dir}{leftgene} = join("-",$leftexons[0],$leftexons[-1]);
    }
    if ($rightexons[0] eq  $rightexons[-1]) {
	$exonnuminfo{$dir}{rightgene} = $rightexons[0];
    }else {
	$exonnuminfo{$dir}{rightgene} = join("-",$rightexons[0],$rightexons[-1]);
    }
}

open OAS, ">$opt{prefix}\.translocations.answer.txt" or die $!;

print OAS join("\t","FusionName","LeftGene","LeftBreakpoint","LeftGeneExons","LeftStrand",
	       "RightGene","RightBreakpoint","RightGeneExons","RightStrand",
	       "RNAReads","DNAReads","FusionType","Annot",'Filter','ChrType','ChrDistance'),"\n";

my $sname = $opt{prefix};

open FUSION, "<$opt{fusion}" or die $!;
my $header = <FUSION>;
chomp($header);
$header =~ s/^#//;
my @hline = split(/\t/,$header);
while (my $line = <FUSION>) {
  chomp($line);
  my @row = split(/\t/,$line);
  my %hash;
  foreach my $i (0..$#row) {
    $hash{$hline[$i]} = $row[$i];
  }
  my @filter;
  my ($left_chr,$left_pos,$left_strand) = split(/:/,$hash{LeftBreakpoint});
  my ($right_chr,$right_pos,$right_strand) = split(/:/,$hash{RightBreakpoint});
  $hash{LeftBreakpoint} = join(":",$left_chr,$left_pos);
  $hash{RightBreakpoint} = join(":",$right_chr,$right_pos);
  $hash{LeftStrand} = $left_strand;
  $hash{RightStrand} = $right_strand;
  $hash{LeftGene} = (split(/\^/,$hash{LeftGene}))[0];
  $hash{RightGene} = (split(/\^/,$hash{RightGene}))[0];
  unless ($keep{$hash{LeftGene}} || $keep{$hash{RightGene}}) {
      push @filter, 'OutsideGeneList';
  }
  $hash{SumRNAReads} += $hash{JunctionReadCount}+$hash{SpanningFragCount};
  my $fname = join("--",$hash{LeftGene},$hash{RightGene});
  my $fname2 = join("--",sort {$a cmp $b} $hash{LeftGene},$hash{RightGene});
  my $key = join("_",$hash{LeftGene},$left_pos,$hash{RightGene},$right_pos);
  my ($leftexon,$rightexon);
  if ($exonnuminfo{$key}) {
      $leftexon = $exonnuminfo{$key}{leftgene};
      $rightexon = $exonnuminfo{$key}{rightgene};
  }
  $hash{PROT_FUSION_TYPE} = 'in-frame' if ($hash{PROT_FUSION_TYPE} eq 'INFRAME');
  my ($dna_support,$rna_support)=("no") x 2;
  $hash{annots} =~ s/CHROMOSOMAL\[/CHROMOSOMAL /;
  $hash{annots} =~ s/\]|\[|\"//g;
  @annots = split(/,/,$hash{annots});
  my ($chrom_type_dist) = grep(/CHROMOSOMAL/,@annots);
  my ($chrtype,$chrdist) = split(/ /,$chrom_type_dist);
  @annot2 = grep(!/CHROMOSOMAL/, @annots);
  $fusion_annot = '';
  if (scalar(@annot2) > 0) {
      $fusion_annot = join(",",@annot2);
  }
  if ($known{$fname2} || $fusion_annot =~ m/CCLE|Cosmic|FA_CancerSupp|Klijn_CellLines|Mitelman|YOSHIHARA_TCGA|chimer/i) {
      push @filter, 'LowReadCt' if ($hash{SumRNAReads} < 3);
  }else {
      push @filter, 'LowReadCt' if ($hash{SumRNAReads} < 5);
  }
  $rna_support = "yes";
  if ($left_chr eq $right_chr) {
      $diff = abs($right_pos-$left_pos);
      push @filter, 'ReadThrough' if ($diff < 200000);
  }
  my $qc ='PASS';
  if (scalar(@filter) > 0) {
      $qc = join(";","FailedQC",@filter);
  }
  
  print OAS join("\t",$fname,$hash{LeftGene},$hash{LeftBreakpoint},$leftexon,$hash{LeftStrand},
		 $hash{RightGene},$hash{RightBreakpoint},$rightexon,$hash{RightStrand},
		 $hash{SumRNAReads},0,lc($hash{PROT_FUSION_TYPE}),$fusion_annot,$qc,$chrtype,$chrdist),"\n";
}


close OAS;
