use strict;

my $listfile = shift;
my $poswfile = shift;
my $burnin = shift;
my $outfile = shift;

open (LISTFILE, $listfile) or die "input error for gene list: $listfile\n";
open (POSWFILE, $poswfile) or die "input error for .posw file: $poswfile\n";
open (OUTFILE, '>'.$outfile) or die "output error: $outfile\n";

my @genename;
my %gene2index;

my $Ngene = <LISTFILE>;
my $ngene = 0;
foreach my $line (<LISTFILE>)	{
	chomp $line;
	$genename[$ngene] = $line;
	$ngene++;
	$gene2index{$line} = $ngene;
}
if ($ngene != $Ngene)	{
	die "error: non matching number of genes: $ngene vs $Ngene\n";
}
	
my @pp;
for (my $i=0; $i<$Ngene; $i++)	{
	@pp[$i] = 0;
}

my $n = 0;
foreach my $line (<POSWFILE>)	{
	$n++;
	if ($n > $burnin)	{
	chomp $line;
	my @a = split('\t',$line);
	my $ngene = (@a);
	if ($ngene != $Ngene)	{
		die "error when reading .posw file: non matching number of genes: $ngene vs $Ngene\n";
	}
	for (my $i=0; $i<$Ngene; $i++)	{
		if ($a[$i] > 0)	{
			$pp[$i] ++;
		}
	}
	}
}


my %gene2pp;
for (my $i=0; $i<$Ngene; $i++)	{
	@pp[$i] /= $n-$burnin;
	$gene2pp{$genename[$i]} = $pp[$i];
}

my $totpp = 0;
my $count = 0;
foreach my $gene (sort {$gene2pp{$b} <=> $gene2pp{$a}} keys %gene2pp)	{
	$count++;
	$totpp += $gene2pp{$gene};
	my $fdr = int(100 * $totpp / $count);
	print OUTFILE "$fdr\t$gene2index{$gene}\t$gene\t$gene2pp{$gene}\n";
}
	


