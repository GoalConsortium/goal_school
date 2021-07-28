#!/usr/bin/perl 
#parse_pindel

my $pair_id = shift @ARGV;
my $vcf = shift @ARGV;

open SI, ">$pair_id.indel.vcf" or die $!;
open SV, ">$pair_id.sv.vcf" or die $!;
open DUP, ">$pair_id.dup.vcf" or die $!;

open VCF, "gunzip -c $vcf|" or die $!;
while (my $line = <VCF>) {
  chomp($line);
  if ($line =~ m/#/) {
    print SI $line,"\n";
    print SV $line,"\n";
    print DUP $line,"\n";
    if ($line =~ m/#CHROM/) {
      my @header = split(/\t/,$line);
      ($chrom, $pos,$id,$ref,$alt,$score,
       $filter,$info,$format,@subjacc) = split(/\t/, $line);
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
  my @deschead = split(/:/,$format);
  my $newformat = 'GT:DP:AD:AO:RO';
  my @newgts = ();
  my $missingGT = 0;
  my @allele;
 FG:foreach my $i (0..$#gts) {
    my $sid = $subjacc[$i];
    my @gtinfo = split(/:/,$gts[$i]);
    my %gtdata;
    if ($allele_info eq '.') {
      push @newgts, '.:.:.:.:.';
      $missingGT ++;
      next FG;
    }
    foreach my $i (0..$#deschead) {
      $gtdata{$deschead[$i]} = $gtinfo[$i];
    }
    if ($gtdata{AD}){
      ($gtdata{RO},@alts) = split(/,/,$gtdata{AD});
      $gtdata{AO} = join(",",@alts);
      $gtdata{DP} = $gtdata{RO};
      foreach (@alts) {
	$gtdata{DP} += $_;
      }
    } elsif (exists $gtdata{NR} && exists $gtdata{NV}) {
      $gtdata{DP} = $gtdata{NR}; 	
      $gtdata{AO} = $gtdata{NV};
      $gtdata{RO} = $gtdata{DP} - $gtdata{AO};
    } elsif (exists $gtdata{AO} && exists $gtdata{RO}) {
      $gtdata{AD} = join(',',$gtdata{RO},$gtdata{AO});
      $gtdata{DP} = $gtdata{RO};
      foreach (split(',',$gtdata{AO})) {
	$gtdata{DP} += $_;
      }
    }
    if (($gtdata{DP} && $gtdata{DP} < 20) || ()) {
      $missingGT ++;
    } if ($gtdata{DP} == 0 || $gtdata{GT} eq './.') {
      push @newgts, '.:.:.:.:.';
      $missingGT ++;
      next FG;
    }
    push @allele, sprintf("%.4f",$gtdata{AO}/$gtdata{DP});
    push @newgts, join(":",$gtdata{GT},$gtdata{DP},$gtdata{AD},$gtdata{AO},$gtdata{RO});
  }
  next if ($missingGT == scalar(@gts));
  
  if ($hash{SVTYPE} eq 'DUP:TANDEM') {
    print DUP join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$annot,$newformat,@newgts),"\n";
  }elsif ($hash{SVTYPE} eq 'DEL' || $hash{SVTYPE} eq 'INS') {
    if (abs($hash{SVLEN}) < 30) {
      print SI join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$annot,$newformat,@newgts),"\n";
    }else {
      $newalt = "<".$hash{SVTYPE}.">";
      print SV join("\t",$chrom,$pos,$id,'N',$newalt,$score,$filter,$annot,$newformat,@newgts),"\n";
    }
  }else {
    $newalt = "<".$hash{SVTYPE}.">";
    print SV join("\t",$chrom,$pos,$id,'N',$newalt,$score,$filter,$annot,$newformat,@newgts),"\n";
  }
}
close VCF;
