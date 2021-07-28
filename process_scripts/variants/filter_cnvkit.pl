#!/usr/bin/perl -w
#parse_cnvkit_table.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'input|s=s','help|h');

my $file = $opt{input};
my $sname = (split(/\./,(split(/\//,$file))[-1]))[0];
my $prefix = (split(/\./,$file))[0];

my %cyto;
open CYTO, "<$prefix\.cytoband.bed" or die $!;
while (my $line = <CYTO>) {
    chomp($line);
    my ($chrom,$start,$end,$band) = split(/\t/,$line);
    my $key = $chrom.":".$start."-".$end;
    $band =~ m/^(\w)(.+)/;
    next unless ($1 && $2);
    push @{$cyto{$key}{$1}},$2;
}

open OUT, ">$prefix\.cnv.answer.txt" or die $!;
open CNSO, ">$prefix\.answerplot.cns" or die $!;

print CNSO join("\t","Chromosome","Start","End","Log2","CN"),"\n";
print OUT join("\t","Gene","Chromosome","Start","End","Abberation Type","CN","Score","CytoBand"),"\n";

open CNR, "<$prefix\.cnr" or die $!;
open CNRO, ">$prefix\.answerplot.cnr" or die $!;
print CNRO join("\t","Gene","Chromosome","Start","End","Log2","Depth","Weight"),"\n";
my $header = <CNR>;
chomp($header);
my @colnames = split(/\t/,$header);

while (my $line = <CNR>) {
    chomp($line);
    my @row = split(/\t/,$line);
    my %hash = ();
    foreach my $j (0..$#row) {
	$hash{$colnames[$j]} = $row[$j];
    }
    my $key = $hash{chromosome}.":".$hash{start}."-".$hash{end};
    my $geneids = $hash{gene};
    my %genes;
    if ($geneids =~ m/ensembl_gn/g) {
	my @ids = split(/;|,/,$geneids);
	foreach my $gid (@ids) {
	    my ($key,$value) = split(/=/,$gid);
	    if ($key eq 'ensembl_gn' || $key eq 'identifier') {
		$genes{$value} = 1;
	    }
	}
    } else {
	my @ids = split(/,/,$geneids);
	foreach my $gid (@ids) {
	    next if ($gid =~ /^rs\d+$|^SNP:rs\d+$|^-$|Fusion|MSI:/);
	    my ($gene,@other) = split(/:/,$gid);
	    my ($gname,@loc) = split(/_/,$gene);
	    $genes{$gname} = 1;
	}
    }
    foreach $gene (keys %genes) {
	print CNRO join("\t",$gene,$hash{chromosome},$hash{start},$hash{end},
			$hash{log2},$hash{depth},$hash{weight}),"\n"; 
    }
}

open IN, "<$file" or die $!;
$header = <IN>;
chomp($header);
@colnames = split(/\t/,$header);

while (my $line = <IN>) {
    chomp($line);
    my @row = split(/\t/,$line);
    my %hash = ();
    foreach my $j (0..$#row) {
	$hash{$colnames[$j]} = $row[$j];
    }
    my $key = $hash{chromosome}.":".$hash{start}."-".$hash{end};
    my $geneids = $hash{gene};
    my %genes;
    if ($geneids =~ m/ensembl_gn/g) {
	my @ids = split(/;|,/,$geneids);
	foreach my $gid (@ids) {
	    my ($key,$value) = split(/=/,$gid);
	    if ($key eq 'ensembl_gn' || $key eq 'identifier') {
		$genes{$value} = 1;
	    }
	}
    } else {
	my @ids = split(/,/,$geneids);
	foreach my $gid (@ids) {
	    next if ($gid =~ /^rs\d+$|^SNP:rs\d+$|^-$|Fusion|MSI:/);
	    my ($gene,@other) = split(/:/,$gid);
	    my ($gname,@loc) = split(/_/,$gene);
	    $genes{$gname} = 1;
	}
    }
    my $len = sprintf("%.1f",($hash{end}-$hash{start})/1000);
    print CNSO join("\t",$hash{chromosome},$hash{start},$hash{end},
		    $hash{log2},$hash{cn}),"\n";
    
    if (exists $hash{cn1} && exists $hash{cn2}) {
      next if (($hash{cn} == 2 && $hash{cn1} ne ''  && $hash{cn2} ne '' && $hash{cn1} > 0 && $hash{cn2} > 0) || scalar(keys %genes) < 1);
      next if ($hash{cn} == 2 && $hash{cn1} eq ''  && $hash{cn2} eq '');
    }else {
      next if ($hash{cn} == 2 || scalar(keys %genes) < 1);
    }
    my $abtype = 'cnLOH';
    if ($hash{cn} < 2) {
      if ($hash{cn} <  1 ) {
	$abtype = 'homozygous deletion';
      }else {
	$abtype = 'hemizygous deletion';
      }
    } elsif ($hash{cn} > 2) {
      $abtype = 'gain' ;
      if ($hash{cn} > 6) {
	$abtype = 'amplification';
      }
      if (exists $hash{cn1} && exists $hash{cn2}) {
	if ($hash{cn1} ne '' && $hash{cn2} ne '' && ($hash{cn1} == 0 || $hash{cn2} == 0)) {
	  $abtype.= ' LOH';
	}
      }
    } else {
      $abtype = 'cnLOH';
    }
    foreach $gene (keys %genes) {
      my @cytoband;
      if ($cyto{$key}{'p'}) {
	@nums = sort {$b <=> $a} @{$cyto{$key}{'p'}};
	push @cytoband, 'p'.$nums[0],'p'.$nums[-1];
      } if ($cyto{$key}{'q'}) {
	@nums = sort {$a <=> $b} @{$cyto{$key}{'q'}};
	push @cytoband, 'q'.$nums[0],'q'.$nums[-1];
      } 
      if ($cytoband[0] eq $cytoband[-1]) {
	$cband = $cytoband[0];
      }else {
	$cband = join("-",$cytoband[0],$cytoband[-1]);
      }
      print OUT join("\t",$gene,$hash{chromosome},$hash{start},$hash{end},
		     $abtype,$hash{cn},$hash{weight},$cband),"\n";
    }
  }
close IN;
