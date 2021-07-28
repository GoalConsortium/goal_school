#!/usr/bin/perl -w
#parse_cnvkit_table.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'skip|s=s','help|h','datadir|r=s','known|k=s','prefix|p=s');


my %keep;
open OM, "<$opt{datadir}/cancer.genelist.txt" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $keep{$line} = 1;
}

my %trxkeep;
open TRX, "<$opt{datadir}/primary_transcripts.txt" or die $!;
while (my $line = <TRX>) {
    chomp($line);
    my ($sym,$gene,$trxid,@other) = split(/\t/,$line);
    $trxkeep{$trxid} = 1;
}

my %knownjunc;
open KNOWN, "<$opt{known}" or die $!;
while (my $line = <KNOWN>) {
  chomp($line);
  my ($chrom,$start,$end,$skipinfo,$chr,$a,$b,$knowninfo,$len) = split(/\t/,$line);    
  my ($stand,$readcount,$known_junction,$exons_skipped) = split(/:/,$knowninfo);
  my $key = join("-",$chrom,$start,$end);
  push @{$knownjunc{$key}}, [$chr,$a,$b,$readcount];
}

my %skip;
open IN, "<$opt{skip}" or die $!;
while (my $line = <IN>) {
  chomp($line);
  my ($chrom,$start,$end,$info,$chr,$a,$b,$gene,$len) = split(/\t/,$line);
  next if ($chr eq '.');
  my ($gname,$readcount,$known_junction,$exons_skipped) = split(/:/,$info);
  my ($sym,$trxid,$exonnum,$strand) = split(/\|/,$gene);
  if ($keep{$sym} && $trxkeep{$trxid}) {
    my $key = join("-",$chrom,$start,$end);
    $skip{$key}{strand} = $strand;
    $skip{$key}{readct} = $readcount;
    $skip{$key}{gene} = $sym;	
    push @{$skip{$key}{trxid}{$trxid}}, $exonnum
  }
}

open OUT, ">$opt{prefix}.exonskip.answer.txt" or die $!;
print OUT join("\t","Sample","Gene","Chromosome","Start","End","Abberation Type","ExonSkipReadct","KnownJunctionReadCt","FractionSkippedFirstJunction","Transcript"),"\n";

foreach my $loci (keys %skip) {
  my $strand = $skip{$loci}{strand};
  my @knownjuncs = @{$knownjunc{$loci}};
  if ($strand eq '+') {
    @knownjuncs = sort {$a->[1] <=> $b->[1]} @knownjuncs;
  }else {    
    @knownjuncs = sort {$b->[2] <=> $a->[2]} @knownjuncs;
  }
  $wtct = $knownjuncs[0][3];
  unless ($wtct) {
    $wtct = 0;
  }
  my @trxs;
  foreach my $trxid (keys %{$skip{$loci}{trxid}}) {
    my @exonnums = sort {$a <=> $b} @{$skip{$loci}{trxid}{$trxid}};
    my $trxname;
    if (scalar(@exonnums) > 1) {
      $trxname = $trxid.':'.$exonnums[0].'-'.$exonnums[-1];
    }else {
      $trxname = $trxid.':'.$exonnums[0];
    }
    print OUT join("\t",$opt{prefix},$skip{$loci}{gene},split(/-/,$loci),'Exon Skipping',
		   $skip{$loci}{readct},$wtct,
		   sprintf("%.4f",$skip{$loci}{readct}/($wtct+$skip{$loci}{readct})),
		   $trxname),"\n";
  }
}
