#!/usr/bin/perl -w
#parse_gencode.pl

my $keep = shift @ARGV;
open KEEP, "<$keep" or die $!;
while (my $line = <KEEP>) {
    chomp($line);
    my ($sym) = split(/\t/,$line);
    $keep{$sym} = 1;
}

open OUT, ">$keep\.bed" or die $!;
open GCODE, "</project/shared/bicf_workflow_ref/GRCh38/gencode.gtf" or die $!;
while (my $line = <GCODE>) {
    chomp($line);
    next if ($line =~ m/^#/);
    my ($chrom,$source,$feature,$start,$end,$filter,$phase,$frame,$info) = 
	split(/\t/,$line);
    next unless ($feature eq 'gene');
    $info =~ s/\"//g;
    my %hash;
    foreach $a (split(/;\s*/,$info)) {
	my ($key,$val) = split(/ /,$a);
	$hash{$key} = $val;
    }
    $hash{gene_id} =~ s/\.\d+//;
    next unless ($keep{$hash{gene_name}});
    print OUT join("\t",$chrom,$start,$end,join("|",$hash{gene_name},$hash{gene_id})),"\n";
}

