#!/usr/bin/perl -w
#integrate_datasets.pl

#module load vcftools/0.1.14 samtools/1.6 bedtools/2.26.0 
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'td|d=s','indel|i=s','sv|s=s','tumor|t=s');

my @files = grep(/vcf.gz/,values %opt);
foreach $file (@files) {
  chomp($file);
  open IN, "gunzip -c $file |" or die $!;
  my $outfile = $file;
  $outfile =~ s/vcf.*/pass.vcf/;
  my @gtheader;
  open OUT, ">$outfile" or die $!;
 W1:while (my $line = <IN>) {
    chomp($line);
    my $format = 'GT:MAF:AO';
    if ($line =~ m/^#/) {
      if ($line =~ m/^#CHROM/) {
	print OUT "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">\n";
	print OUT "##FORMAT=<ID=MAF,Number=1,Type=Integer,Description=\"Mutation Allele Frequency\">\n";
	print OUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	my @header = split(/\t/,$line);
	($chrom, $pos,$id,$ref,$alt,$score,
	 $filter,$info) = split(/\t/, $line);
	print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,
		       $filter,$info,'FORMAT',$opt{tumor}),"\n";
      }else {
	print OUT $line,"\n";
      }
      next;
    }
    my ($chrom, $pos,$id,$ref,$alt,$score,
	$filter,$annot) = split(/\t/, $line);
    next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
    my %hash = ();
    foreach $a (split(/;/,$annot)) {
      my ($key,$val) = split(/=/,$a);
      $hash{$key} = $val unless ($hash{$key});
    }
    next unless ($hash{ANN});
    next unless ($hash{ANN} =~ m/HIGH|MODERATE|LOW/);
    my %gtinfo;
    foreach (split(/,/,$hash{DP2})) {
      $gtinfo{AO} += $_;
    }
    $gtinfo{MAF} = $hash{VAF};
    $hash{AF} = $hash{VAF};
    next if $gtinfo{AO} < 10;
    next if ($hash{VAF} < 0.01);
    $gt = '0/0';
    $gt = '0/1' if ($hash{VAF} > 0.3);
    $gt = '1/1' if ($hash{VAF} > 0.8);
    my $gtinfo = join(":",$gt,$gtinfo{MAF},$gtinfo{AO});
    my @nannot;
    foreach $info (sort {$a cmp $b} keys %hash) {
      if (defined $hash{$info}) {
	push @nannot, $info."=".$hash{$info};
      }else {
	push @nannot, $info;
      }
    }
    my $newannot = join(";",@nannot);
    print OUT join("\t",$chrom, $pos,$id,$ref,$alt,$score,$filter,$newannot,
		   $format,$gtinfo),"\n";
  }
  close IN;
}
