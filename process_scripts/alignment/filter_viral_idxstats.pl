#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;
use Getopt::Long;
#use Cwd;

#my $cwd = getcwd();

my $pair_id = "";
GetOptions ('pairid|p=s' => \$pair_id);
if($pair_id eq ""){$pair_id = "all";}

my @allFiles = @ARGV;
chomp(@allFiles);

my %data;
my @samples;
open OUT, ">$pair_id\.viral.seqstats.txt" or die $!;

my %virus;
$virus{'NC_001538.1'}="BK polyomavirus";
$virus{'NC_007605.1'}="Human gammaherpesvirus 4";
$virus{'NC_009334.1'}="Human herpesvirus 4";
$virus{'NC_001488.1'}="Human T-lymphotropic virus 2";
$virus{'D90400.1'}="Human papillomavirus type 58";
$virus{'M14119.1'}="Human papillomavirus type 11 (HPV-11)";
$virus{'J04353.1'}="Human papillomavirus type 31 (HPV-31)";
$virus{'M12732.1'}="Human papillomavirus type 33";
$virus{'M74117.1'}="Human papillomavirus type 35";
$virus{'M62849.1'}="Human papillomavirus ORFs";
$virus{'X74481.1'}="Human papillomavirus type 52";
$virus{'X77858.1'}="Human papilloma virus type 59";
$virus{'NC_001355.1'}="Human papillomavirus type 6b";
$virus{'NC_001357.1'}="Human papillomavirus - 18";
$virus{'NC_001436.1'}="Human T-lymphotropic virus 1";
$virus{'NC_001699.1'}="JC polyomavirus";
$virus{'EF177177.1'}="Human papillomavirus type 56 clone Qv26342";
$virus{'NC_009333.1'}="Human herpesvirus 8";
$virus{'EF202167.1'}="Human papillomavirus type 45 isolate Qv31748";
$virus{'NC_014407.1'}="Human polyomavirus 7";
$virus{'HM355825.1'}="Merkel cell polyomavirus isolate MCVw156";
$virus{'NC_001526.4'}="Human papillomavirus type 16";


print OUT join("\t","SampleName","VirusAcc","VirusName","Mapped","Unmapped"),"\n";
foreach my $file_idx(@allFiles){
		print "File Included: ".$file_idx."\n";
		#my @filePath = split("/", $file_idx);
		my $fileName = $file_idx;
		$fileName =~ s/\.idxstats\.txt//g;
		push @samples, $fileName;

		open INFILE, "<$file_idx" or die $!;

		foreach my $line(<INFILE>){
				chomp $line;
				my @data_array = split("\t",$line);
				if($data_array[2]!=0 and $data_array[0] ne '*'){
					$data{$data_array[0]}{$fileName} = $data_array[2]."\t".$data_array[3];
					print OUT join("\t",$fileName,$data_array[0],$virus{$data_array[0]},$data_array[2],$data_array[3]),"\n";
				}
		}
		close INFILE;
}
