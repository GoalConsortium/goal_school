#!/usr/bin/perl -w
#integrate_datasets.pl

#module load vcftools/0.1.14 samtools/1.6 bedtools/2.26.0 
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'in|i=s','pid|p=s','tumor|t=s');

open VCFOUT, ">$opt{pid}\.vcf" or die $!;
open DELFUS, ">$opt{pid}\.genefusion.txt" or die $!;

open IN, "gunzip -c $opt{in} |" or die $!;
my @gtheader;
W1:while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#/) {
    if ($line =~ m/^#CHROM/) {
      print VCFOUT "##INFO=<ID=AF,Number=A,Type=Integer,Description=\"Alternate allele observation frequency\">\n";
      my @header = split(/\t/,$line);
      my ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@sids) = split(/\t/, $line);
      @gtheader = @sids;
      print DELFUS join("\t",'CHROM','POS','CHR2','END','Effect',
			'Gene','BioType','Filter','Format',@gtheader),"\n";
      unless ($opt{tumor}) {
	if (grep(/T_DNA/,@gtheader)) {
	  my @tsamps = grep(/T_DNA/,@gtheader);
	  $opt{tumor} = $tsamps[0];
	}else {
	  $opt{tumor} = $gtheader[0];
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
  next unless ($hash{ANN});
  next unless ($hash{ANN} =~ m/HIGH|MODERATE|LOW/);
  unless ($hash{SVTYPE}) {
      $hash{SVTYPE} = $hash{VTYPE};
  }
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
  if ($filter ne 'PASS') {
    $filter = 'FailedQC;'.$filter;
  }
  my ($svid,$svpair) = split(/:/,$id);
  if ($svid =~ m/(\w+_\d+)(\w)$/) {
      ($svid,$svpair) = ($1,$2);
      $svpair =~ tr/oh/12/;
  }
  $svpair = 1 unless ($svpair);

  my $newline = join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$newannot,$format,@gts);
  my ($allele,$effect,$impact,$gene,$geneid,$feature,
      $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
      $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$keeptrx);
  unless ($hash{CHR2}) {
      $hash{CHR2} = $chrom;
  }
  my $fusionline = join("\t",$chrom,$pos,$hash{CHR2},$hash{END},$effect,$gene,$biotype,$filter,$format,@gts);
  $svpairs{$svid}{$svpair} = {vcfline=>$line,fusionline=>$fusionline,
			      filter=>$filter,gene=>$keepforvcf,
			      svtype=>$hash{SVTYPE},
			      alt=>$alt,span=>$hash{SPAN}};
}
close IN;

foreach my $id (keys %svpairs) {
  my $alt1 = $svpairs{$id}{1}{alt};
  my $svtype = $svpairs{$id}{1}{svtype};
  if (scalar(keys %{$svpairs{$id}}) > 1) {
    my $alt2 = $svpairs{$id}{2}{alt};
    if ($alt1 =~ m/^\w\[/ && $alt2 =~ m/^\]/) {
      $svtype = 'DEL';
    }elsif ($alt2 =~ m/^\w\[/ && $alt1 =~ m/^\]/) {
      $svtype = 'INS';
    }
  }
  if (($svtype eq 'DEL' && $svpairs{$id}{1}{span} && $svpairs{$id}{1}{span} > 9999)
      || $svpairs{$id}{1}{fusionline} =~ m/fusion|transcript_ablation/
      || $svpairs{$id}{1}{gene} =~ m/&/) {
    print DELFUS $svpairs{$id}{1}{fusionline},"\n";
  } elsif ($svtype eq 'DUP' || $svtype eq 'INS' || ($svtype eq 'DEL' && $svpairs{$id}{1}{gene} !~ m/&/ && $svpairs{$id}{1}{span} < 9999)) {
    print VCFOUT $svpairs{$id}{1}{vcfline},"\n"
  }
}
