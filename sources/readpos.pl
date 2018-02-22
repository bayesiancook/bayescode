use strict;

my $infile = shift;
my $listfile = shift;
my $burnin = shift;

open (LIST0, $listfile) or die "input error for original gene list\n";
open (LISTFILE, $infile.'.genelist') or die "input error for gene list\n";
open (POSWFILE, $infile.'.posw') or die "input error for .posw file\n";
open (POSOMFILE, $infile.'.posom') or die "input error for .posom file\n";

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
my @meanw;
my @meanom;
for (my $i=0; $i<$Ngene; $i++)	{
	@pp[$i] = 0;
	@meanw[$i] = 0;
	@meanom[$i] = 0;
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
            $meanw[$i] += $a[$i];
        }
	}
}

$n = 0;
foreach my $line (<POSOMFILE>)	{
	$n++;
	if ($n > $burnin)	{
        chomp $line;
        my @a = split('\t',$line);
        my $ngene = (@a);
        if ($ngene != $Ngene)	{
            die "error when reading .posw file: non matching number of genes: $ngene vs $Ngene\n";
        }
        for (my $i=0; $i<$Ngene; $i++)	{
            $meanom[$i] += $a[$i];
        }
	}
}

my %gene2pp;
my %gene2meanom;
my %gene2meanw;
for (my $i=0; $i<$Ngene; $i++)	{
	$pp[$i] /= $n-$burnin;
	$meanw[$i] /= $n-$burnin;
	$meanom[$i] /= $n-$burnin;

	$gene2pp{$genename[$i]} = $pp[$i];
	$gene2meanw{$genename[$i]} = $meanw[$i];
	$gene2meanom{$genename[$i]} = $meanom[$i];
}

open (SOUTFILE, '>'.$infile.'.sortedpp') or die "output error: $infile\n";
my $totpp = 0;
my $count = 0;
foreach my $gene (sort {$gene2pp{$b} <=> $gene2pp{$a}} keys %gene2pp)	{
	$count++;
	$totpp += $gene2pp{$gene};
	my $fdr = int(100 * $totpp / $count);
	print SOUTFILE "$fdr\t$gene2pp{$gene}\t$gene2meanw{$gene}\t$gene2meanom{$gene}\t$gene\n";
}
	
open (PPOUTFILE, '>'.$infile.'.pp') or die "output error: $infile\n";
for (my $g=0; $g<$Ngene; $g++)  {
    my $gene = $genename0[$g];
    print PPOUTFILE "$gene\t$gene2pp{$gene}\t$gene2meanw{$gene}\t$gene2meanom{$gene}\n";
}


