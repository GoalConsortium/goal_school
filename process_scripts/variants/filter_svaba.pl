#!/usr/bin/perl -w
#integrate_datasets.pl

#module load vcftools/0.1.14 samtools/1.6 bedtools/2.26.0 
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'in|i=s','pid|p=s','tumor|t=s','sv|s=s',
			  'method|a=s');

open VCFOUT, ">$opt{pid}\.$opt{method}.vcf" or die $!;
open DELFUS, ">$opt{pid}\.$opt{method}.potentialfusion.txt" or die $!;

open IN, "gunzip -c $opt{in} |" or die $!;

W1:while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#/) {
    if ($line =~ m/^#CHROM/) {
      print VCFOUT "##INFO=<ID=AF,Number=A,Type=Integer,Description=\"Alternate allele observation frequency\">\n";
      my @header = split(/\t/,$line);
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@gtheader) = split(/\t/, $line);
      unless ($opt{tumor}) {
	if (grep(/T_DNA/,@gtheader)) {
	  my @tsamps = grep(/T_DNA/,@gtheader);
	  $opt{tumor} = $tsamps[0];
	}else {
	  $opt{tumor} = $gtheader[0];
	}
      }
            unless ($opt{normal}) {
	if (grep(/N_DNA/,@gtheader)) {
	  my @tsamps = grep(/N_DNA/,@gtheader);
	  $opt{normal} = $tsamps[0];
	}
      }
     }
    print VCFOUT $line,"\n";
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
  next unless ($chrom =~ m/^chr\d+$/ || $chrom eq 'chrX' || $chrom eq 'chrY');
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val unless ($hash{$key});
  }
  if (length($alt) > length($ref)) {
    $hash{SVTYPE} = "INS";
    $hash{END} = $pos;
  }elsif (length($alt) < length($ref)) {
    my $diff = substr($ref, length($alt));
    $hash{SVTYPE} = "DEL";
      $hash{END} = $pos + length($diff);
  }
  next unless ($hash{ANN});
  next unless ($hash{ANN} =~ m/HIGH|MODERATE|LOW/);
  my %gtinfo = ();
  my @deschead = split(/:/,$format);
 F1:foreach my $k (0..$#gtheader) {
    my $subjid = $gtheader[$k];
    my $allele_info = $gts[$k];
    my @ainfo = split(/:/, $allele_info);
    my @mutallfreq = ();
    foreach my $k (0..$#ainfo) {
      $gtinfo{$subjid}{$deschead[$k]} = $ainfo[$k];
    }
    $gtinfo{$subjid}{DP} = (split(/,/,$gtinfo{$subjid}{DP}))[0] if ($gtinfo{$subjid}{DP});
    next F1 unless ($gtinfo{$subjid}{DP} && $gtinfo{$subjid}{DP} ne '.' && $gtinfo{$subjid}{DP} >= 1);
    my @altct = split(/,/,$gtinfo{$subjid}{AD});
    my $refct = shift @altct;
    @altct2 = split(/,/,$gtinfo{$subjid}{AO});
    if (scalar(@altct) ne scalar(@altct2)) {
      warn "Inconsistent Allele counts @ $chrom,$pos,$alt,$gtinfo{$subjid}{AD},$gtinfo{$subjid}{AO}\n";
    }
    my $total = $refct;
    foreach  my $act (@altct) {
      next if ($act eq '.');
      $total += $act;
      push @mutallfreq, sprintf("%.4f",$act/$gtinfo{$subjid}{DP});
    }
    $gtinfo{$subjid}{MAF} = \@mutallfreq;
  }
  next unless ($gtinfo{$opt{tumor}}{DP} && $gtinfo{$opt{tumor}}{DP} ne '.' && $gtinfo{$opt{tumor}}{DP} >= 20);
  unless ($gtinfo{$opt{tumor}}{AO} =~ m/\d+/ && $gtinfo{$opt{tumor}}{AD} =~ m/,/) {
    warn "Missing Alt:$line\n";
  }
  @tumormaf = @{$gtinfo{$opt{tumor}}{MAF}};
  @tumoraltct = split(/,/,$gtinfo{$opt{tumor}}{AO});
  next if ($tumoraltct[0] eq '.');
  $hash{AF} = join(",",@tumormaf);
  next if ($tumoraltct[0] < 5);
  my $keepforvcf = 0;
  my $keeptrx;
 F1:foreach $trx (split(/,/,$hash{ANN})) {
    my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	$cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
    next unless ($impact =~ m/HIGH|MODERATE|LOW/ || $effect =~ /splice/i);
    #next if($effect eq 'sequence_feature');
    $keeptrx = $trx;
    $keepforvcf = $gene;
    last F1;
  }
  next unless $keepforvcf;
  my @nannot;
  foreach $info (sort {$a cmp $b} keys %hash) {
    if (defined $hash{$info}) {
      push @nannot, $info."=".$hash{$info};
    }else {
      push @nannot, $info;
    }
  }
  my $newannot = join(";",@nannot);
  if ($hash{SVTYPE} eq 'INS' || ($hash{SVTYPE} eq 'DEL' && $keepforvcf !~ m/&/)) {
    if ($filter =~ m/LOWMAPQ|LowQual/i) {
      $filter = 'FailedQC;'.$filter;
    }
    print VCFOUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$newannot,
		      $format,@gts),"\n";
  } elsif ($hash{SVTYPE} eq 'DEL' && $hash{SPAN} && $hash{SPAN} > 9999) {
    my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	$cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$keeptrx);
    print DELFUS join("\t",$chrom,$pos,$hash{CHR2},$hash{END},$effect,$gene,$biotype,$filter,$format,@gts),"\n";
  }
}
close IN;

my %svpairs;

open IN, "gunzip -c $opt{sv} |" or die $!;

W1:while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#/) {
    if ($line =~ m/^#CHROM/) {
      my @header = split(/\t/,$line);
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@gtheader) = split(/\t/, $line);
      unless ($opt{tumor}) {
	if (grep(/T_DNA/,@gtheader)) {
	  my @tsamps = grep(/T_DNA/,@gtheader);
	  $opt{tumor} = $tsamps[0];
	}else {
	  $opt{tumor} = $gtheader[0];
	}
      }
    }
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val unless ($hash{$key});
  }
  if ($alt =~ m/CHR(\w+):(\d+)/i) {
    $hash{CHR2} = 'chr'.$1;
    $hash{END} = $2;
  }

  next unless ($hash{ANN});
  next unless ($hash{ANN} =~ m/HIGH|MODERATE|LOW/);
  my %gtinfo = ();
  my @deschead = split(/:/,$format);
 F1:foreach my $k (0..$#gtheader) {
    my $subjid = $gtheader[$k];
    my $allele_info = $gts[$k];
    my @ainfo = split(/:/, $allele_info);
    my @mutallfreq = ();
    foreach my $k (0..$#ainfo) {
      $gtinfo{$subjid}{$deschead[$k]} = $ainfo[$k];
    }
    $gtinfo{$subjid}{DP} = (split(/,/,$gtinfo{$subjid}{DP}))[0] if ($gtinfo{$subjid}{DP});
    next F1 unless ($gtinfo{$subjid}{DP} && $gtinfo{$subjid}{DP} ne '.' && $gtinfo{$subjid}{DP} >= 1);
    my @altct = split(/,/,$gtinfo{$subjid}{AD});
    my $refct = shift @altct;
    @altct2 = split(/,/,$gtinfo{$subjid}{AO});
    if (scalar(@altct) ne scalar(@altct2)) {
      warn "Inconsistent Allele counts @ $chrom,$pos,$alt,$gtinfo{$subjid}{AD},$gtinfo{$subjid}{AO}\n";
    }
    my $total = $refct;
    foreach  my $act (@altct) {
      next if ($act eq '.');
      $total += $act;
      push @mutallfreq, sprintf("%.4f",$act/$gtinfo{$subjid}{DP});
    }
    $gtinfo{$subjid}{MAF} = \@mutallfreq;
  }
  next unless ($gtinfo{$opt{tumor}}{DP} && $gtinfo{$opt{tumor}}{DP} ne '.' && $gtinfo{$opt{tumor}}{DP} >= 20);
  unless ($gtinfo{$opt{tumor}}{AO} =~ m/\d+/ && $gtinfo{$opt{tumor}}{AD} =~ m/,/) {
    warn "Missing Alt:$line\n";
  }
  @tumormaf = @{$gtinfo{$opt{tumor}}{MAF}};
  @tumoraltct = split(/,/,$gtinfo{$opt{tumor}}{AO});
  next if ($tumoraltct[0] eq '.');
  $hash{AF} = join(",",@tumormaf);
  next if ($tumoraltct[0] < 5);
  #next if ($tumormaf[0] < 0.01);
  my $keepforvcf = 0;
  my $keeptrx;
 F1:foreach $trx (split(/,/,$hash{ANN})) {
    my ($allele,$effect,$impact,$gene,$geneid,$feature,
	$featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
	$cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
    next unless ($impact =~ m/HIGH|MODERATE/ || $effect =~ /splice/i);
    #next if($effect eq 'sequence_feature');
    $keeptrx = $trx;
    $keepforvcf = $gene;
    last F1;
  }
  next unless $keepforvcf;
  my @nannot;
  foreach $info (sort {$a cmp $b} keys %hash) {
    if (defined $hash{$info}) {
      push @nannot, $info."=".$hash{$info};
    }else {
      push @nannot, $info;
    }
  }
  my $newannot = join(";",@nannot);
  my ($svid,$svpair) = split(/:/,$id);
  my $newline = join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$newannot,$format,@gts);
  my ($allele,$effect,$impact,$gene,$geneid,$feature,
      $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
      $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$keeptrx);
  my $fusionline = join("\t",$chrom,$pos,$hash{CHR2},$hash{END},$effect,$gene,$biotype,$filter,$format,@gts);
  $svpairs{$svid}{$svpair} = {vcfline=>$line,fusionline=>$fusionline,filter=>$filter,
			      gene=>$keepforvcf,alt=>$alt,span=>$hash{SPAN}};
}
close IN;


foreach my $id (keys %svpairs) {
  my $alt1 = $svpairs{$id}{1}{alt};
  my $alt2 = $svpairs{$id}{2}{alt};
  my $svtype;
  if ($alt1 =~ m/^\w\[/ && $alt2 =~ m/^\]/) {
    $svtype = 'DEL';
  }elsif ($alt2 =~ m/^\w\[/ && $alt1 =~ m/^\]/) {
    $svtype = 'INS';
  }else {
    $svtype = 'UNK';
  }
  if (($svtype eq 'DEL' && $svpairs{$id}{1}{span} && $svpairs{$id}{1}{span} > 9999) ||($svpairs{$id}{1}{fusionline} =~ m/fusion/)) {
      print DELFUS $svpairs{$id}{1}{fusionline},"\n";
  } elsif ($svtype eq 'INS' || ($svtype eq 'DEL' && $svpairs{$id}{1}{gene} !~ m/&/ && $svpairs{$id}{1}{span} < 9999)) {
      if ($svpairs{$id}{1}{filter} =~ m/LOWMAPQ|LowQual/i) {
	  my @vline = split(/\t/,$svpairs{$id}{1}{vcfline});
	  $vline[6] = 'FailedQC'.$vline[6];
	  $svpairs{$id}{1}{vcfline} = join("\t",@vline);
      }
      print VCFOUT $svpairs{$id}{1}{vcfline},"\n"
  }
}
