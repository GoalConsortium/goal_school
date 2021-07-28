#!/usr/bin/perl 
#migrate_db.pl

my $pair_id = shift @ARGV;
my $vcf = shift @ARGV;
my $outfile = $pair_id.".uniform.vcf";
open OUT, ">$outfile" or die $!;
open VCF, "gunzip -c $vcf|" or die $!;
my @headerline;
while (my $line = <VCF>) {
  chomp($line);
  if ($line =~ m/#/) {
    next if ($line =~ m/FORMAT=<ID=AO/ || $line =~ m/FORMAT=<ID=AD/ || $line =~ m/FORMAT=<ID=RO/ || $line =~ m/FORMAT=<ID=DP/);
    if ($line =~ m/#CHROM/) {
      print OUT "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">\n";
      print OUT "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">\n";
      print OUT "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
      print OUT "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n";
      print OUT "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n";
      print OUT "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n";
      print OUT "##INFO=<ID=VTYPE,Number=1,Type=String,Description=\"Variant Type\">\n";
      print OUT "##INFO=<ID=SPAN,Number=1,Type=Integer,Description=\"Variant Size\">\n";
      @headerline = split(/\t/, $line);
      my ($c, $p,$i,$r,$a,$s,$f,$an,$fo,@snames) = @headerline;
      foreach my $j (0..$#snames) {
	$snames[$j] =~ s/\[|\]|\.consensus|\.final//g;
      }
      print OUT join("\t",$c, $p,$i,$r,$a,$s,$f,$an,$fo,@snames),"\n";
    } else {
      print OUT $line,"\n";
    }
    next;
  }
  my @line = split(/\t/, $line);
  my $numfields = scalar(@line);
  my ($chrom, $pos,$id,$ref,$alt,$score,$filter,$annot,$format) = @line[0..8];
  if ($numfields != scalar(@headerline)) {
    my @snames = @headerline[9..$#headerline];
    my $numsamples = scalar(@snames);
    my $idx = $numfields-$numsamples;
    @gts = @line[$idx..$#line];
  } else {
    @gts = @line[9..$#line];
  }
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val;
  }
  $refnt = $ref;
  $altnt = $alt;
  if (lc($alt) =~ m/^(\w+)[\]|\[](\w+):(\d+)[\]|\[]/) {
    ($altnt,$chr2,$end) = ($1,$2,$3);
    $hash{END} = $end;
    $hash{CHR2} = $chr2;
  }elsif (lc($alt) =~ m/^[\]|\[](\w+):(\d+)[\]|\[](\w+)/) {
    ($chr2,$end,$altnt) = ($1,$2,$3);
    $hash{END} = $end;
    $hash{CHR2} = $chr2;
  }
  if ($altnt =~ m/^\w+$/ && $refnt  =~ m/^\w+$/) {
    if (length($altnt) > length($refnt)) {
      $hash{VTYPE} = "INS";
      $hash{END} = $pos unless $hash{END};
    }elsif (length($altnt) < length($refnt)) {
      my $diff = substr($ref, length($alt));
      $hash{VTYPE} = "DEL";
      $hash{END} = $pos + length($diff) unless $hash{END};
    }
  }
  if ($hash{END}) {
    $hash{SPAN} = abs($hash{END} - $pos);
  }
  if ($hash{INSERTION}) {
    $hash{VTYPE} = "INS";
  }
  if ($hash{SVTYPE} && $hash{VTYPE} && $hash{SVTYPE} eq 'BND') {
    $hash{SVTYPE} = $hash{VTYPE};
  }
  my @deschead = split(/:/,$format);
  my $newformat = 'GT:DP:AD:AO:RO';
  my @newgts = ();
  my $missingGT = 0;
  #next if (scalar(@gts) != scalar(@snames));
 FG:foreach my $allele_info (@gts) {
    my @gtinfo = split(/:/,$allele_info);
    my %gtdata;
    if ($allele_info eq '.') {
      push @newgts, '.:.:.:.:.';
      $missingGT ++;
      next FG;
    }
    foreach my $i (0..$#deschead) {
      $gtdata{$deschead[$i]} = $gtinfo[$i];
    }
    if ($gtdata{AD} =~ m/\d+,\d+/){
      ($gtdata{RO},@alts) = split(/,/,$gtdata{AD});
      $gtdata{AO} = join(",",@alts);
      $gtdata{DP} = $gtdata{RO};
      foreach (@alts) {
	$gtdata{DP} += $_;
      }
    } elsif ($gtdata{AD} =~ m/^\d+$/){
      $gtdata{AO} = $gtdata{AD};
      $gtdata{RO} = $gtdata{DP} - $gtdata{AO};
      if ($gtdata{RO} < 0) {
	$gtdata{DP} +=  $gtdata{AO};
	$gtdata{RO} = $gtdata{DP} -  $gtdata{AO};
      }
      $gtdata{AD} = join(',',$gtdata{RO},$gtdata{AO});
    } elsif (exists $gtdata{DV} && exists $gtdata{RV}) {
      $gtdata{AO} = $gtdata{DV} + $gtdata{RV};
      $gtdata{RO} = $gtdata{DR} + $gtdata{RR};
      $gtdata{AD} = join(',',$gtdata{RO},$gtdata{AO});
      $gtdata{DP} = $gtdata{RO}+$gtdata{AO};
    } elsif (exists $gtdata{DR} && exists $gtdata{SR}){
      $gtdata{AO} = $gtdata{AD};
      $gtdata{DP} = $gtdata{AO} unless $gtdata{DP};
      if  ($gtdata{DP} > $gtdata{AD}) {
	$gtdata{RO} = $gtdata{DP} - $gtdata{AD};
      } else {
	$gtdata{RO} = 0;
      }
      $gtdata{AD} = join(',',$gtdata{RO},$gtdata{AO});
    } elsif (exists $gtdata{NR} && exists $gtdata{NV}) {
      $gtdata{DP} = $gtdata{NR}; 	
      $gtdata{AO} = $gtdata{NV};
      $gtdata{RO} = $gtdata{DP} - $gtdata{AO};
      $gtdata{AD} = join(',',$gtdata{RO},$gtdata{AO});
    } elsif (exists $gtdata{AO} && exists $gtdata{RO}) {
      $gtdata{AD} = join(',',$gtdata{RO},$gtdata{AO});
      $gtdata{DP} = $gtdata{RO};
      foreach (split(',',$gtdata{AO})) {
	$gtdata{DP} += $_;
      }
    } elsif ($gtdata{TIR}) {
      $gtdata{GT} = '0/0';
      $gtdata{AO} = (split(/,/,$gtdata{TIR}))[0];
      $gtdata{RO} = $gtdata{DP} - $gtdata{AO};
      $gtdata{AD} = join(',',$gtdata{RO},$gtdata{AO});
    } elsif ($gtdata{VF}) {
      $gtdata{GT} = '0/0';
      $gtdata{AO} = $gtdata{VF};
      $gtdata{DP} = $gtdata{AO} unless $gtdata{DP};
      if  ($gtdata{DP} > $gtdata{AO}) {
	$gtdata{RO} = $gtdata{DP} - $gtdata{AD};
      } else {
	$gtdata{RO} = 0;
      }
      $gtdata{AD} = join(',',$gtdata{RO},$gtdata{AO});
    } elsif ($gtdata{$ref."U"} && $gtdata{$alt."U"}) {
      $gtdata{GT} = '0/0';
      $gtdata{AO} = (split(/,/,$gtdata{$alt."U"}))[0];
      $gtdata{RO} = (split(/,/,$gtdata{$ref."U"}))[0];
      $gtdata{AD} = join(',',$gtdata{RO},$gtdata{AO});
    }
    if ($gtdata{DP} && $gtdata{DP} < 5) {
      $missingGT ++;
    }
    if ($gtdata{DP} == 0 || $gtdata{GT} eq './.') {
      push @newgts, '.:.:.:.:.';
      $missingGT ++;
      next FG;
    }
    push @newgts, join(":",$gtdata{GT},$gtdata{DP},$gtdata{AD},$gtdata{AO},$gtdata{RO});
  }
  next if ($missingGT == scalar(@gts));
  my @nannot;
  foreach $info (sort {$a cmp $b} keys %hash) {
    if (defined $hash{$info}) {
      push @nannot, $info."=".$hash{$info};
    }else {
      push @nannot, $info;
    }
  }
  my $newannot = join(";",@nannot);
  print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$newannot,$newformat,@newgts),"\n";
}
close VCF;
