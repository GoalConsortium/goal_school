#!/usr/bin/perl -w
#integrate_datasets.pl

my $bamreadct = shift @ARGV;
open NRC, "<$bamreadct" or die $!;
open OUT, ">$bamreadct\.cttable.txt" or die $!;
while (my $line = <NRC>) {
  chomp($line);
  my ($chrom,$pos,$ref,$depth,@reads) = split(/\t/,$line);
  next unless ($depth > 10);
  $chrom = 'chr'.$chrom if ($chrom !~ m/^chr/);
  my $ro;
  my %hash;
  foreach my $rct (@reads) {
    my ($base,$basect,@otherstats) = split(/:/,$rct);
    if ($ref eq $base) {
      $ro = $basect;
    }else {
      if ($base =~ m/\+|\-/) {
	$base =~ s/\+/$ref/;
	#$base =~ s/\-/$ref/;
      }
      $hash{$base} = $basect if ($basect);
    }
  }
  my @basecalls;
  foreach (keys %hash) {
    push @basecalls, join(":",$_,$hash{$_});
  }
  print OUT join("\t",$chrom,$pos,$depth,$ref,$ro,join(";",@basecalls)),"\n";
}
close NRC;
