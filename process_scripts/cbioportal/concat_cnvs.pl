#!/usr/bin/perl -w
#concat_cnvs.pl

my @answercnv = @ARGV;

my %cts;
my %sample;
foreach $file (@answercnv) {
    open IN, "<$file" or die $!;
    my $fname = (split(/\//,$file))[-1];
    my ($sample,@ext) = split(/\./,$fname);
    $sample{$sample} = 1;
    my $header = <IN>;
    chomp($header);
    my @cols = split(/\t/, $header);
    while (my $line = <IN>) {
	chomp($line);
	my %hash;
	my @row = split(/\t/,$line);
	foreach my $i (0..$#row) {
	    $hash{$cols[$i]} = $row[$i];
	}
	next if ($hash{Chromosome} =~ m/X|Y/);
	$ct = $hash{CN};
	
	my $discreet = 0;
	if ($ct == 1) {
	    $discreet = -1;
	}elsif ($ct == 0) {
	    $discreet = -2;
	}elsif ($ct == 3) {
	    $discreet = 1;
	}elsif ($ct > 3) {
	    $discreet = 2;
	}elsif ($ct eq 'NA') {
	    $ct = '';
	}
	$gene=$hash{Gene};
	$discreet{$gene}{$sample} = $discreet;
	$cn{$gene}{$sample} = $hash{CN};
    }
}
my @samples = sort {$a cmp $b} keys %sample;
open DIS, ">discreet.cna.txt" or die $!;
open CN, ">copynumber.cna.txt" or die $!;
print DIS join("\t",'Hugo_Symbol',@samples),"\n";
print CN join("\t",'Hugo_Symbol',@samples),"\n";
foreach my $gene (keys %discreet) {
    my @line1;
    my @line2;
    foreach my $sid (@samples) {
	$discreet{$gene}{$sid} = 0 unless (exists $discreet{$gene}{$sid});
	$cn{$gene}{$sid} = 2 unless (exists $cn{$gene}{$sid});
	push @line1, $discreet{$gene}{$sid};
	push @line2, $cn{$gene}{$sid};
    }
    print DIS join("\t",$gene,@line1),"\n";
    print CN join("\t",$gene,@line2),"\n";
}
