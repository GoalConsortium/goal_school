#!/usr/bin/perl -w
#concat_cnvs.pl

my @discreet = `ls *cnv_discreet.txt`;
my %cts;
my %sample;
foreach $file (@discreet) {
    open IN, "<$file" or die $!;
    my ($sample,@ext) = split(/\./,$file);
    $sample{$sample} = 1;
    while (my $line = <IN>) {
	chomp($line);
	my ($chr,$s,$e,$ct,$gene) = split(/\t/,$line);
	my $discreet = 0;
	if ($ct eq 'NA') {
	    $ct = '';
	}elsif ($ct == 1) {
	    $discreet = -1;
	}elsif ($ct == 0) {
	    $discreet = -2;
	}elsif ($ct == 3) {
	    $discreet = 1;
	}elsif ($ct > 3) {
	    $discreet = 2;
	}
	$cts{$gene}{$sample} = $discreet;
    }
}
my @samples = sort {$a cmp $b} keys %sample;
open OUT, ">discreet.cna.txt" or die $!;
print OUT join("\t",'Hugo_Symbol',@samples),"\n";
foreach my $gene (keys %cts) {
    my @line;
    foreach my $sid (@samples) {
	$cts{$gene}{$sid} = 0 unless (exists $cts{$gene}{$sid});
	push @line, $cts{$gene}{$sid};
    }
    print OUT join("\t",$gene,@line),"\n";
}

my @continuous = `ls *cnv_continuous.txt`;
my %cts;
my %sample;
foreach $file (@continuous) {
    open IN, "<$file" or die $!;
    my ($sample,@ext) = split(/\./,$file);
    $sample{$sample} = 1;
    while (my $line = <IN>) {
	chomp($line);
	my ($chr,$s,$e,$ct,$gene) = split(/\t/,$line);
	$cts{$gene}{$sample} = $ct;
    }
}
my @samples = sort {$a cmp $b} keys %sample;
open OUT, ">continuous.cna.txt" or die $!;
print OUT join("\t",'Hugo_Symbol',@samples),"\n";
foreach my $gene (keys %cts) {
    my @line;
    foreach my $sid (@samples) {
	$cts{$gene}{$sid} = 2 unless ($cts{$gene}{$sid});
	push @line, $cts{$gene}{$sid};
    }
    print OUT join("\t",$gene,@line),"\n";
}
