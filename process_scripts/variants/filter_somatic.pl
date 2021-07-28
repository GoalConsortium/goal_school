#!/usr/bin/perl -w
#integrate_datasets.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'prefix|p=s','normal|n=s','dnavcf|f=s',
			  'help|h');


open OUT, ">$opt{prefix}\.filt.vcf" or die $!;
my @sampids;
open IN, "gunzip -c $opt{dnavcf} |" or die $!;
W1:while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    my ($chrom, $pos,$id,$ref,$alt,$score,
     $filter,$info,$format,@gtheader) = split(/\t/, $line);
    @sampids = @gtheader;
    print OUT join("\t",$line),"\n";
    next;
  } elsif ($line =~ m/^#/) {
    print OUT $line,"\n";
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  my %fail;
  next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val unless ($hash{$key});
  }
  my @exacaf;
  my $exacaf;
  if ($hash{AF_POPMAX}) {
    foreach (split(/,/,$hash{AF_POPMAX})) {
      push @exacaf, $_ if ($_ ne '.');
    }
    @exacaf = sort {$b <=> $a} @exacaf;
    $exacaf = $exacaf[0] if ($exacaf[0]);
  } if ($hash{dbNSFP_ExAC_Adj_AF}) {
    foreach (split(/,/,$hash{dbNSFP_ExAC_Adj_AF})) {
      push @exacaf, $_ if ($_ ne '.');
    }
    @exacaf = sort {$b <=> $a} @exacaf;
    if ($exacaf[0]) {
      if ($exacaf && $exacaf[0] < $exacaf ) {
	$exacaf[0] = $exacaf;
      }else {
	$exacaf = $exacaf[0] if ($exacaf[0]);
      }
    } 
  }elsif ($hash{AC_POPMAX} && $hash{AN_POPMAX}) {
    my @exacs = split(/,/,$hash{AC_POPMAX});
    my $ac = 0;
    foreach $val (@exacs) {
      $ac += $val if ($val =~ m/^\d+$/);
    }
    my @exans = split(/,/,$hash{AN_POPMAX});
    my $an = 0;
    foreach $val (@exans) {
      $an += $val if ($val =~ m/^\d+$/);
    }
    $exacaf = sprintf("%.4f",$ac/$an) if ($ac > 0 && $an > 10);
  }
  next if ($exacaf && $exacaf > 0.05);
  $fail{'COMMON'} = 1 if ($exacaf && $exacaf > 0.01);
  $fail{'StrandBias'} = 1 if (($hash{FS} && $hash{FS} > 60) || $filter =~ m/strandBias/i || $hash{strandBias} || (($hash{SAP} && $hash{SAP} > 20) && ((exists $hash{SAF} && $hash{SAF}< 1 && $hash{SRF} > 10) || (exists $hash{SAR} && $hash{SAR}< 1 && $hash{SRR} > 10))));
  my $cosmicsubj = 0;
  if ($hash{CNT}) {
    my @cosmicct = split(/,/,$hash{CNT}); 
    foreach $val (@cosmicct) {
      $cosmicsubj += $val if ($val =~ m/^\d+$/);
    }
  }
  my %gtinfo = ();
  my @tumoraltct;
  my @tumormaf;
  my @tumordepth;
  my @deschead = split(/:/,$format);
  $hash{SS} = 5  unless ($hash{SS});
 F1:foreach my $k (0..$#sampids) {
    my $subjid = $sampids[$k];
    my $allele_info = $gts[$k];
    my @ainfo = split(/:/, $allele_info);
    my @mutallfreq = ();
    foreach my $k (0..$#ainfo) {
      $gtinfo{$deschead[$k]} = $ainfo[$k];
      $hash{$deschead[$k]} = $ainfo[$k] unless ($opt{normal} && $subjid eq $opt{normal});
    }
    $gtinfo{DP} = (split(/,/,$gtinfo{DP}))[0] if ($gtinfo{DP});
    next F1 unless ($gtinfo{DP} && $gtinfo{DP} ne '.' && $gtinfo{DP} >= 1);
    
    my @altct = split(/,/,$gtinfo{AO});
    
    foreach  my $act (@altct) {
      next if ($act eq '.');
      push @mutallfreq, sprintf("%.4f",$act/$gtinfo{DP});
    }
    if ($opt{normal} && $subjid eq $opt{normal}) {
      $hash{NormalAF} = $mutallfreq[0];
      $hash{NormalDP} = $gtinfo{DP};
      if ($mutallfreq[0] >= 0.25) {
	$hash{SS} = 1;
      }elsif ($mutallfreq[0] >= 0.05) {
	$hash{'HighFreqNormalAF'} = 1;
      }else {
	$hash{SS} = 2;
      }
    } else {
      push @tumoraltct, @altct;
      push @tumordepth, $gtinfo{DP};
      push @tumormaf, @mutallfreq;
    }
  }
  @tumoraltct = sort {$b <=> $a} @tumoraltct;
  @tumormaf = sort {$b <=> $a} @tumormaf;
  @tumordepth = sort {$a <=> $b} @tumordepth;
  next unless ($tumordepth[0] && $tumordepth[0] ne '.' && $tumordepth[0] >= 20);
  if ($hash{NormalAF} && $hash{NormalAF} < 0.05 &&  $hash{NormalAF}*5 > $tumormaf[0]) {
    $hash{'HighFreqNormalAF'} = 1;
    $hash{SS} = 5;
  }
  if (exists $hash{INDEL}) {
    $hash{TYPE} = 'indel';
  }
  $hash{TYPE} = 'ambi' unless ($hash{"TYPE"});
  next if ($tumoraltct[0] eq '.');
  $hash{AF} = join(",",@tumormaf);
  my @callers;
  if ($hash{CallSet} && $hash{CallSet} =~ m/\|/) {
    my @callinfo ;
    @callinfo = split(/,/, $hash{CallSet}) if ($hash{CallSet});
    foreach $cinfo (@callinfo) {
      my ($caller, $alt, @samafinfo) = split(/\|/,$cinfo);
      push @callers, $caller;
    }
    $hash{CallSet} = join(",",@callinfo);
  }
  if (($id =~ m/COS/ && $cosmicsubj >= 5) || $hash{OncoKBHotspot}) {
    $fail{'LowAltCt'} = 1 if ($tumoraltct[0] < 3);
    $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.05);
    $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.1 && $hash{TYPE} ne 'snp');
  }else {
    $fail{'OneCaller'} = 1 if (scalar(@callers) < 2);
    $fail{'LowAltCt'} = 1 if ($tumoraltct[0] < 8);
    $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.05);
    $fail{'LowMAF'} = 1 if ($tumormaf[0] < 0.10 && $hash{TYPE} ne 'snp');
  }
  if ($hash{RepeatType} && $hash{RepeatType} =~ m/Simple_repeat/ && $tumormaf[0] < 0.15) {
    $fail{'InRepeat'} = 1;
  }
  delete $hash{SOMATIC};
  my $keepforvcf;
  if ($hash{ANN}) {
      foreach $trx (split(/,/,$hash{ANN})) {
	  my ($allele,$effect,$impact,$gene,$geneid,$feature,
	      $featureid,$biotype,$rank,$codon,$aa,$cdna_pos,$len_cdna,
	      $aapos,$distance,$err) = split(/\|/,$trx);
	  next unless ($impact =~ m/HIGH|MODERATE/ || $effect =~ /splice/i);
	  next if($effect eq 'sequence_feature');
	  $keepforvcf = $gene;
      }
  }
  unless ($keepforvcf) {
    $fail{'NonCoding'} = 1;
  }
  my @fail = sort {$a cmp $b} keys %fail;
  if (scalar(@fail) == 0) {
    $filter = 'PASS';
  }else {
    $filter = join(";", 'FailedQC',@fail);
  }
  my @nannot;
  foreach $info (sort {$a cmp $b} keys %hash) {
    if (defined $hash{$info}) {
      push @nannot, $info."=".$hash{$info};
    }else {
      push @nannot, $info;
    }
  }

  my $newannot = join(";",@nannot);
  print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,$filter,$newannot,$format,@gts),"\n";
}
close IN;
