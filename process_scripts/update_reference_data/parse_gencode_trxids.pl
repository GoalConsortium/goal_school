#!/usr/bin/perl -w
#parse_gencode.pl
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'help|h','protein|p');
my %keep;
open TRX, "<cosmic_transcripts.txt" or die $!;
while (my $line = <TRX>) {
    chomp($line);
    my ($geneid,$enst) = split(/\t/,$line);
    $keep{$enst} = 1;
}

open EBED, ">temp.bed" or die $!;
my $gtf = shift @ARGV;
open GCODE, "$gtf" or die $!;
while (my $line = <GCODE>) {
    chomp($line);
    next if ($line =~ m/^#/);
    my ($chrom,$source,$feature,$start,$end,$phase,$strand,$frame,$info) = 
	split(/\t/,$line);
    $info =~ s/\"//g;
    my %hash;
    foreach $a (split(/;\s*/,$info)) {
	my ($key,$val) = split(/ /,$a);
	$hash{$key} = $val;
    }
    $hash{gene_id} =~ s/\.\d+//;
    next unless $hash{transcript_id};
    next if ($opt{protein} && $hash{gene_type} ne 'protein_coding');
    next unless $keep{(split(/\./,$hash{transcript_id}))[0]};
    if ($feature eq 'CDS') {
	$hash{transcript_id} =~ s/\.\d+//;
	print EBED join("\t",$chrom,$start,$end,join("|",$hash{gene_name},$hash{transcript_id},$hash{exon_number},$strand)),"\n";
    }
}
close EBED;
system("sort -V -k 1,1 -k 2,2n temp.bed  |bedtools merge -o distinct -c 4 -i stdin > cosmictrx.gencode.exon.bed");
