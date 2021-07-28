#!/usr/bin/perl 
#migrate_db.pl

my $pair_id = shift @ARGV;
my $vcf = shift @ARGV;
my $outfile = $pair_id.".vcf";
open OUT, ">$outfile" or die $!;
open VCF, "$vcf" or die $!;
while (my $line = <VCF>) {
    chomp($line);
    if ($line =~ m/#/) {
	if ($line =~ m/#CHROM/) {
	    print OUT "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">\n";
	    print OUT "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">\n";
	    print OUT "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
	    print OUT "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n";
	}
	unless ($line =~ m/FORMAT=<ID=AO/ || $line =~ m/FORMAT=<ID=AD/ || $line =~ m/FORMAT=<ID=RO/ || $line =~ m/FORMAT=<ID=DP/) {
	    print OUT $line,"\n";
	}
	next;
    }
    my ($chrom, $pos,$id,$ref,$alt,$score,
	$filter,$annot,$format,@gts) = split(/\t/, $line);
    next if ($chrom =~ m/alt/);
    my %hash = ();
    foreach $a (split(/;/,$annot)) {
	my ($key,$val) = split(/=/,$a);
	$hash{$key} = $val;
	}
    my @deschead = split(/:/,$format);
    my $newformat = 'GT:DP:AD:AO:RO';
    my @newgts = ();
    my $missingGT = 0;
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
      if ($alt eq '.' || $alt eq '<NON_REF>') {
	  $gtdata{AO} = 0;
	  $gtdata{RO} = $gtdata{DP};
	  $gtdata{AD} = join(",", $gtdata{RO},$gtdata{AO});
      } elsif ($gtdata{AD}){
	  ($gtdata{RO},@alts) = split(/,/,$gtdata{AD});
	  $gtdata{AO} = join(",",@alts);
	  $gtdata{DP} = $gtdata{RO};
	  foreach (@alts) {
	      $gtdata{DP} += $_;
	  }
	  
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
    if ($hash{END}) {
	foreach $i ($pos..$hash{END}) {
	    print OUT join("\t",$chrom,$i,$id,$ref,'.',$score,$filter,$annot,$newformat,@newgts),"\n";
	}
    }else {
	print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$annot,$newformat,@newgts),"\n";
    }
}
close VCF;
