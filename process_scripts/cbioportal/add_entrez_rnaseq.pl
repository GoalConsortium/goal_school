#!/usr/bin/perl -w
#patient_sample_uuid.pl

my $studyid = shift @ARGV;
open ENT_ENS, "</project/shared/bicf_workflow_ref/human/gene_info.human.txt" or die $!;
my %entrez;
my $ent_header = <ENT_ENS>;
while (my $line = <ENT_ENS>){
  chomp $line;
  my @row = split(/\t/, $line);
  $entrez{$row[2]}=$row[1];
}
close ENT_ENS;

my %rnaids;

open INF, "<countTable.fpkm.txt" or die $!;
open OUTF, ">expression.txt" or die $!;
my $inheader = <INF>;
chomp($inheader);
my @incol = split(/\t/,$inheader);
my @newcol = ();
foreach my $cname (@incol) {
    next if ($cname =~ m/ENSEMBL|TYPE/);
    if ($cname =~ m/SYMBOL/) {
	push @newcol, ('Hugo_Symbol','Entrez_Gene_Id');
    }else {
	push @newcol, $cname;
	$rnaids{$cname} = 1;
    }
}
print OUTF join("\t",@newcol),"\n";
while (my $line = <INF>) {
    chomp($line);
    my ($ens,$sym,$type,@nums) = split(/\t/,$line);
    next if ($type ne 'protein_coding');
    my ($min,$max) = min(@nums);
    next unless ($max > 1);
    next unless ($entrez{$sym});
    print OUTF join("\t",$sym,$entrez{$sym},@nums),"\n";
}

close INF;
close OUTF;

open RNAI, ">case_lists/rnaseq.txt" or die $!;
print RNAI join("\n","cancer_study_identifier: $studyid","stable_id: $studyid\_rna_seq_mrna",
		"case_list_name: mRNA","case_list_description: RNAseq Samples",
		"case_list_ids:".join("\t",keys %rnaids)),"\n";
close RNAI;

sub sum {
        my @data = @_;
        my $total = 0;
        foreach (@data) {
	    $total += $_;
        }
        return $total;
}
sub min {
        my @data = @_;
	@data = sort {$a <=> $b} @data;
        return ($data[0],$data[-1]);
}
