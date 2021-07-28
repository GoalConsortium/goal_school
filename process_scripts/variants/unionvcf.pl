#!/usr/bin/perl 
#migrate_db.pl

my $headerfile = shift @ARGV;
my @vcffiles = @ARGV;
my @callers = ('fb','mutect','gatk','platypus','strelka2','shimmer','pindel','svaba');
%algos = map {$_=>1} @callers;

open HEADER, "<$headerfile" or die $!;
open OUT, ">int.vcf" or die $!;
while (my $line = <HEADER>) {
  chomp($line);
  print OUT $line,"\n";;
}
close HEADER;
my @sampleorder;

my %headerlines;
foreach $vcf (@vcffiles) {
  my @filename = (split(/\./,$vcf));
  my $caller;
  foreach $fio (@filename) {
    $caller = $fio if ($algos{$fio});
  }
  open VCF, "gunzip -c $vcf|" or die $!;
  my @sampleids;
  while (my $line = <VCF>) {
    chomp($line);
    if ($line =~ m/#/) {
      if ($line =~ m/#CHROM/) {
	($chromhd, $posd,$idhd,$refhd,$althd,$scorehd,
	 $filterhd,$annothd,$formathd,@sampleids) = split(/\t/, $line);
	foreach $j (0..$#sampleids) {
	  $sampleids[$j] = (split(/\./,$sampleids[$j]))[0];
	}
	unless (@sampleorder) {
	  @sampleorder = @sampleids;
	  print OUT join("\t",$chromhd, $posd,$idhd,$refhd,$althd,$scorehd,
			 $filterhd,$annothd,$formathd,@sampleids),"\n";
	}
	next;
      }
      next;
    }
    my ($chrom, $pos,$id,$ref,$alt,$score,
	$filter,$annot,$format,@gts) = split(/\t/, $line);
    my %hash = ();
    foreach $a (split(/;/,$annot)) {
      my ($key,$val) = split(/=/,$a);
      $hash{$key} = $val;
    }
    #next if (($hash{FS} && $hash{FS} > 60) || $filter =~ m/strandBias/i);
    my @deschead = split(/:/,$format);
    my $newformat = 'GT:DP:AD:AO:RO';
    my %newgts;
    my %afinfo;
    my %gtfilt;
    my $missingGT = 0;
  FG:foreach my $i (0..$#gts) {
      my $allele_info = $gts[$i];
      my $sid = $sampleids[$i];
      my @gtinfo = split(/:/,$allele_info);
      my %gtdata;
      if ($allele_info eq '.') {
	$newgts{$sid} = '.:.:.:.:.';
	$missingGT ++;
	next FG;
      }
      foreach my $i (0..$#deschead) {
	$gtdata{$deschead[$i]} = $gtinfo[$i];
      }
      if ($gtdata{GT} eq './.') {
	$newgts{$sid} = '.:.:.:.:.';
	$missingGT ++;
	next FG;
      }
      if ($gtdata{FT} && $gtdata{FT} =~ m/HighSNVSB/) {
	  $gtfilt{'StrandBias'} = 1;
      }
      if ($gtdata{DP4}) { #varscan uses this
	my ($ref_fwd,$ref_rev,$alt_fwd,$alt_rev) = split(',',$gtdata{DP4});
	$gtdata{AO} = $alt_fwd+$alt_rev;
	$gtdata{RO} = $ref_fwd+$ref_rev;
	$gtdata{DP} = $ref_fwd+$ref_rev+$alt_fwd+$alt_rev;
	$gtdata{AD} = join(",",$gtdata{RO},$gtdata{AO});
      }elsif ($gtdata{AD} && $gtdata{AD} =~ m/,/){
	($gtdata{RO},@alts) = split(/,/,$gtdata{AD});
	$gtdata{AO} = join(",",@alts);
	$gtdata{DP} = $gtdata{RO};
	foreach (@alts) {
	  $gtdata{DP} += $_;
	}
      } elsif (exists $gtdata{NR} && exists $gtdata{NV}) { #platypus uses this
	$gtdata{DP} = (split(/,/,$gtdata{NR}))[0]; 	
	$gtdata{AO} = $gtdata{NV};
	$gtdata{RO} = $gtdata{DP};	
	foreach $altct (split(/,/,$gtdata{NV})) {
	  $gtdata{RO} -= $altct;
	}
	$gtdata{AD} = join(",",$gtdata{RO},$gtdata{AO})
      } elsif (exists $gtdata{AO} && exists $gtdata{RO}) {
	$gtdata{AD} = join(',',$gtdata{RO},$gtdata{AO});
	$gtdata{DP} = $gtdata{RO};
	foreach (split(',',$gtdata{AO})) {
	  $gtdata{DP} += $_;
	}
      }
      my $mafs = '.';
      my @maf = ();
      if ($gtdata{DP} && $gtdata{DP} ne '.' && exists $gtdata{AO} && $gtdata{AO} ne '.') {
	foreach $areadct (split(/,/,$gtdata{AO})) {
	  push @maf, sprintf("%.2f",$areadct/$gtdata{DP});
	}
	$mafs = join(",",@maf);
      }
      if (exists $gtdata{DP} && $gtdata{DP} < 20) {
	$missingGT ++;
      }elsif (exists $gtdata{AO} && $gtdata{AO} < 3) {
	$missingGT ++;
      }
      $afinfo{$sid} = join(":",$gtdata{DP},$mafs);
      $newgts{$sid} =  join(":",$gtdata{GT},$gtdata{DP},$gtdata{AD},$gtdata{AO},$gtdata{RO});
    }
    next if ($missingGT == scalar(@gts));
    my @gtdesc;
    my @newgts;
    foreach $id (@sampleorder) {
      push @gtdesc, join(":",$id,$afinfo{$id});
      push @newgts, $newgts{$id};
    }
    my @filts = split(";",$filter);
    my %filterqc = map {$_ => 1} @filts;
    delete $filterqc{'.'};
    if ($gtfilt{'StrandBias'} || ($hash{FS} && $hash{FS} > 60)) { # ||($hash{SAP} && $hash{SAP} > 20)) {
	$filterqc{strandBias} = 1;
	$annot .= ';strandBias=1';
    }if (scalar(keys %filterqc) > 0) {
	$filter = join(";",keys %filterqc);
    }else {
	$filter = '.';
    }
    $lines{$chrom}{$pos}{$ref}{$alt}{$caller} = [$chrom,$pos,$id,$ref,$alt,$score,$filter,$annot,$newformat,\@newgts,\@gtdesc] unless $lines{$chrom}{$pos}{$ref}{$alt}{$caller};
  }
  close VCF;
}

F1:foreach $chr (sort {$a cmp $b} keys %lines) {
 F2:foreach $pos (sort {$a <=> $b} keys %{$lines{$chr}}) {
  F5:foreach $ref (sort {$a <=> $b} keys %{$lines{$chr}{$pos}}) {
    F4:foreach $alt (sort {$a <=> $b} keys %{$lines{$chr}{$pos}{$ref}}) {
	my @callset;
	my %csets;
	my %depth;
      F3:foreach $caller (sort {$a cmp $b} keys %{$lines{$chr}{$pos}{$ref}{$alt}}) {
	  my ($chrom, $pos,$id,$ref,$alt,$score,$filter,$annot,
	      $format,$gtsref,$gtdescref) = @{$lines{$chr}{$pos}{$ref}{$alt}{$caller}};
	  my @gtdesc = @{$gtdescref};
	  my @gtdesc2;	
	  foreach $gtd (@gtdesc) {
	    my ($id,$dp,$maf) = split(/:/,$gtd);
	    push @gtdesc2, $dp;
	    push @{$csets{$id}}, [$caller,$dp,$maf];
	  }
	  @gtdesc2 = sort {$b <=> $a} @gtdesc2;
	  $depth{$caller} = $gtdesc2[0];
	  push @callset, join("|",$caller,$alt,@gtdesc);
	}
	my $consistent = 1;
	foreach $id (keys %csets) {
	  my @calls = @{$csets{$id}};
	  my @calls = sort {$a[2] <=> $b[2]} @calls;
	  $consistent = 0 if ($calls[0][2] < 0.25 && $calls[-1][2] - $calls[0][2] > 0.10 && $calls[-1][2]/($calls[0][2]+0.001) > 3);
	}
	@callorder = sort {$depth{$b} <=> $depth{$a}} keys %depth;
      F3:foreach $caller (@callers) {
	  if ($lines{$chr}{$pos}{$ref}{$alt}{$caller}) {
	    my ($chrom, $pos,$id,$ref,$alt,$score,$filter,$annot,
		$format,$gtsref,$gtdescref) = @{$lines{$chr}{$pos}{$ref}{$alt}{$caller}};
	    @gts = @{$gtsref};
	    @gtdesc = @{$gtdescref};
	    $annot = $annot.";CallSet=".join(",",@callset);
	    unless ($consistent) {
	      $annot = $annot.";CallSetInconsistent=1";
	    }
	    print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,
			   $filter,$annot,$format,@gts),"\n";
	    last F3;
	  }
	}
      }
    }
  }
}
close OUT;
