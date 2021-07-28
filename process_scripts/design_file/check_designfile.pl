#!/usr/bin/perl -w
#check_designfile.pl

use strict;
use warnings;

my $pe = shift @ARGV;
my $dfile = shift @ARGV;
open OUT, ">design.valid.txt" or die $!;
open DFILE, "<$dfile" or die $!;
my $head = <DFILE>;
chomp($head);
$head =~ s/FullPathTo//g;
my @colnames = split(/\t/,$head);
my %newcols = map {$_=> 1} @colnames;

unless (grep(/FqR1/,@colnames)) {
    die "Missing Sequence Files in Design File: FqR1\n";
}
unless (grep(/SampleID/,@colnames)) {
    die "Missing SampleID in Design File\n";
}

if ($pe eq 'pe') {
    unless (grep(/FqR2/,@colnames)) {
	die "Missing Sequence Files in Design File: FqR2\n";
    }
}else {
    delete $newcols{FqR2};
}

my @cols = sort {$a cmp $b} keys %newcols;
print OUT join("\t",@cols),"\n";
my @grp = ('a','b');
my $lnct = 0;
while (my $line = <DFILE>) {
    chomp($line);
    $line =~ s/FullPathTo//g;
    my @row = split(/\t/,$line);
    my %hash;
    foreach my $i (0..$#row) {
	next unless ($newcols{$colnames[$i]});
	$row[$i] =~ s/-//g unless ($colnames[$i] =~ m/Fq/);
	$hash{$colnames[$i]} = $row[$i];
    }
    if ($hash{SampleID} =~ m/^\d/) {
	$hash{SampleID} =~ s/^/S/;
    }
    $hash{SampleName} = $hash{SampleID} unless ($hash{SampleName});
    $hash{SubjectID} = $hash{SampleID} unless ($hash{SubjectID});
    unless ($hash{SampleGroup}) {
	my $j = $lnct %% 2;
	$hash{SampleGroup} = $grp[$j];
    }
    $lnct ++;
    $hash{SampleGroup} =~ s/_//g;
    next unless ( -f $hash{FqR1});
    if ($hash{FqR1} =~ m/gz/) {
	system(qq{mv $hash{FqR1} $hash{SampleID}.R1.fastq.gz});
    }else {
	system(qq{mv $hash{FqR1} $hash{SampleID}.R1.fastq});
	system(qq{pigz -f $hash{SampleID}.R1.fastq});
    }
    $hash{FqR1} = "$hash{SampleID}.R1.fastq.gz";
    if ($pe eq 'pe') {
	next unless (-f $hash{FqR2});
	if ($hash{FqR2} =~ m/gz/) {
	    system(qq{mv $hash{FqR2} $hash{SampleID}.R2.fastq.gz});
	}else {
	    system(qq{mv $hash{FqR2} $hash{SampleID}.R2.fastq});
	    system(qq{pigz -f $hash{SampleID}.R2.fastq});
	}
	$hash{FqR2} = "$hash{SampleID}.R2.fastq.gz";
    }
    my @line;
    foreach my $f (@cols) {
	push @line, $hash{$f};
    }
    print OUT join("\t",@line),"\n";
    print join(",",$hash{SampleID},"$hash{SampleID}.R1.fastq.gz","$hash{SampleID}.R2.fastq.gz"),"\n";
}
