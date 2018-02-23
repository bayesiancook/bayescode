use strict;

my $infile = shift;
my $listfile = shift;
my $burnin = shift;

open (LIST0, $listfile) or die "input error for original gene list\n";
open (LISTFILE, $infile.'.genelist') or die "input error for gene list\n";
open (OMFILE, $infile.'.geneom') or die "input error for .posom file\n";

my @genename0;
my @genename;
my %gene2size;

my $Ngene0 = <LIST0>;
my $ngene0 = 0;
foreach my $line (<LIST0>)	{
	chomp $line;
	$genename0[$ngene0] = $line;
	$ngene0++;
}
if ($ngene0 != $Ngene0)	{
	die "error: non matching number of genes: $ngene0 vs $Ngene0\n";
}

my $Ngene = <LISTFILE>;
if ($Ngene != $Ngene0)  {
	die "error: non matching number of genes: $Ngene vs $Ngene0\n";
}

my $ngene = 0;
foreach my $line (<LISTFILE>)	{
	chomp $line;
    my @a = split('\t',$line);
    my $n = (@a);
    if ($n != 2)    {
        die "error when parsing genelist file: $line\n";
    }
    my $name = $a[0];
    my $size = $a[1];
    $genename[$ngene] = $name;
    $gene2size{$name} = $size;
	$ngene++;
}
if ($ngene != $Ngene)	{
	die "error: non matching number of genes: $ngene vs $Ngene\n";
}
	
my @pp;
my @meanom;
for (my $i=0; $i<$Ngene; $i++)	{
	@pp[$i] = 0;
	@meanom[$i] = 0;
}

my $n = 0;
foreach my $line (<OMFILE>)	{
	$n++;
	if ($n > $burnin)	{
        chomp $line;
        my @a = split('\t',$line);
        my $ngene = (@a);
        if ($ngene != $Ngene)	{
            die "error when reading .posw file: non matching number of genes: $ngene vs $Ngene\n";
        }
        for (my $i=0; $i<$Ngene; $i++)	{
            if ($a[$i] > 1.0)	{
                $pp[$i] ++;
            }
            $meanom[$i] += $a[$i];
        }
	}
}

my %gene2pp;
my %gene2meanom;
for (my $i=0; $i<$Ngene; $i++)	{
	$pp[$i] /= $n-$burnin;
	$meanom[$i] /= $n-$burnin;

	$gene2pp{$genename[$i]} = $pp[$i];
	$gene2meanom{$genename[$i]} = $meanom[$i];
}

open (SOUTFILE, '>'.$infile.'.sortedpostom') or die "output error: $infile\n";
my $totpp = 0;
my $count = 0;
foreach my $gene (sort {$gene2meanom{$b} <=> $gene2meanom{$a}} keys %gene2meanom)	{
	print SOUTFILE "$gene2meanom{$gene}\t$gene2pp{$gene}\t$gene\n";
}
	
open (PPOUTFILE, '>'.$infile.'.postom') or die "output error: $infile\n";
for (my $g=0; $g<$Ngene; $g++)  {
    my $gene = $genename0[$g];
    print PPOUTFILE "$gene\t$gene2pp{$gene}\t$gene2meanom{$gene}\n";
}


