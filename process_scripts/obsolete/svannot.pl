#!/usr/bin/perl -w
#svanno.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'input|i=s','refdb|r=s','help|h');

my $vcf = $opt{input};
my $outfile = $vcf;
$outfile =~ s/vcf/txt/g;
open OUT, ">$outfile" or die $!;
my $gffile = $vcf;
$gffile =~ s/vcf/genefusion.txt/g;
open GF, ">$gffile" or die $!;

my %lines;
my %eventid;
my $ct = 0;
open IN, "<$vcf" or die $!;
while (my $line = <IN>) {
  chomp($line);
  if ($line =~ m/^#CHROM/) {
    my @header = split(/\t/,$line);
    ($chrom, $pos,$id,$ref,$alt,$score,
     $filter,$info,$format,@subjacc) = split(/\t/, $line);
  }
  next if $line =~ m/#/;
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  $ct ++;
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val unless ($hash{$key});
  }
  if ($id eq 'N') {
      $id = 'NB'.sprintf("%06s",$ct);
  }
  my $evid = (split(/_/,$id))[0];
  $hash{'END'} = $pos+1 unless $hash{'END'};
  my ($allele,$effect,$impact,$gene,$geneid,$feature,
      $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
      $cds_pos,$cds_len,$aapos,$aalen,$distance,$err);
  next unless $hash{ANN};
 F1:foreach $trx (split(/,/,$hash{ANN})) {
     ($allele,$effect,$impact,$gene,$geneid,$feature,
      $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
      $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
     next unless ($impact && $impact =~ m/HIGH|MODERATE/);
     next if ($effect eq 'sequence_feature');
     last F1;
 }
  next unless ($impact && $impact =~ m/HIGH|MODERATE/);
  next unless ($gene);
  next if ($done{$chrom}{$pos});
  $done{$chrom}{$pos} = 1;
  @deschead = split(":",$format);
 F1:foreach $sample (@subjacc) {
     my $allele_info = shift @gts;
     @ainfo = split(/:/, $allele_info);
     my %gtinfo = ();
     foreach $k (0..$#deschead) {
	 $gtinfo{$deschead[$k]} = $ainfo[$k];
     }
     unless ($gtinfo{SU}) {
	 $gtinfo{SU} = 0;
	 $gtinfo{SU} = $gtinfo{RV}+$gtinfo{DV} if ($gtinfo{RV} && $gtinfo{DV});
     }
     next if ($gtinfo{SU} < 10); 
     if (lc($alt) =~ m/chr(\w+):(\d+)/i) {
	 $chr2 = $1;
	 $end = $2;
	 if ($chr2 eq $chrom) {
	     print OUT join("\t",$sample,$gene,$chrom,$pos,$end,$id,$hash{SVTYPE},$effect,
			    $featureid,$codon,$aa,$rank,$gtinfo{SU}),"\n";
	 }elsif ($id =~ m/_\d+/ && $effect =~ m/gene_fusion/) {
	     $gene =~ s/\&/--/;
	     my ($left_gene,$right_gene) = split(/--/,$gene);
	     next unless ($right_gene);
	     my $left_breakpoint = join(":",$chrom,$pos);
	     my $right_breakpoint = join(":",'chr'.$chr2,$end);
	     my ($leftstrand,$rightstrand) = split(//,$hash{STRANDS});
	     print GF join("\t",$sample,$gene,$left_gene,$right_gene, $left_breakpoint,
			    $right_breakpoint,$leftstrand,$rightstrand,'',$gtinfo{SU}),"\n";
	 }else {
	     print OUT join("\t",$sample,$gene,$chrom,$pos,$end,$ id,$hash{SVTYPE},$effect,
			    $featureid,$codon,$aa,$rank,$gtinfo{SU}),"\n";
	 }
     }else {
	 if ($hash{CHR2} && $hash{CHR2} eq $chrom) {
	     $end = $hash{END};
	     print OUT join("\t",$sample,$gene,$chrom,$pos,$end,$id,$hash{SVTYPE},$effect,
			    $featureid,$codon,$aa,$rank,$gtinfo{SU}),"\n";
	 }elsif ($hash{CHR2} && $hash{CHR2} ne $chrom) {
	     $gene =~ s/\&/--/;
	     my ($left_gene,$right_gene) = split(/--/,$gene);
	     my $left_breakpoint = join(":",$chrom,$pos);
	     my $right_breakpoint = join(":",$hash{CHR2},$hash{END});
	     print GF join("\t",$sample,$gene,$left_gene,$right_gene, $left_breakpoint,$right_breakpoint,'','','',$gtinfo{SU}),"\n";
	 }unless ($hash{CHR2}) {
	     $hash{END} = $pos + 1 unless ($hash{END});
	     print OUT join("\t",$sample,$gene,$chrom,$pos,$hash{END},$id,$hash{SVTYPE},$effect,
			    $featureid,$codon,$aa,$rank,$gtinfo{SU}),"\n";
	 }
     }
 }
}

close IN;
