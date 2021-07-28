#!/usr/bin/perl -w
#merge_tables.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;

my %opt = ();
my $results = GetOptions (\%opt,'outdir|o=s','help|h');

$opt{outdir} = './' unless ($opt{outdir});

my @files = @ARGV;
foreach $file (@files) {
    chomp($file);
    open IN, "<$file" or die $!;
    $fname = basename($file);
    my $sample = (split(/\./,$fname))[0];
    $samps{$sample} = 1;
    my $head = <IN>;
    chomp($head);
    my @colnames = split(/\t/,$head);
    while (my $line = <IN>) {
	chomp($line);
	my @row = split(/\t/,$line);
	my $gene = $row[0];
	my $ct = $row[-1];
	$cts{$gene}{$sample} = $ct;
    }
    close IN;
}

open OUT, ">$opt{outdir}\/countFeature.txt" or die $!;
my @samples = sort {$a cmp $b} keys %samps;
print OUT join("\t",'FeatureCtType',@samples),"\n";
foreach $gene (keys %cts) {
    my @line;
    foreach $sample (@samples) {
	unless ($cts{$gene}{$sample}) {
	    $cts{$gene}{$sample} = 0;
	}
	push @line, $cts{$gene}{$sample};
    }
    print OUT join("\t",$gene,@line),"\n";
}
close OUT;
