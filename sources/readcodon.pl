
use strict;

my $basename = shift;
my $burnin = shift;
my $cutoff = shift;
my $outfile = shift;

open (LISTFILE, $basename.'.genelist') or die "input error for gene list\n";

my $Ngene = <LISTFILE>;
chomp $Ngene;

my %gene2nsite;
my %gene2offset;
my @genename;
my $totnsite = 0;

my $gene = 0;
foreach my $line (<LISTFILE>)   {
    chomp $line;
    my @a = split('\t',$line);
    my $n = (@a);
    if ($n != 2)    {
        die "error when reading gene: $line\n";
    }
    my $name = $a[0];
    my $nsite = $a[1];
    $gene2nsite{$name} = $nsite;
    $genename[$gene] = $name;
    $gene2offset{$name} = $totnsite;
    $totnsite += $nsite;
    $gene++;
}
if ($gene != $Ngene) {
    die "error: non matching number of genes in list file\n";
}

open (POSWFILE, $basename.'.posw') or die "input error for pow file\n";

my @pp;
for (my $i=0; $i<$Ngene; $i++)	{
	@pp[$i] = 0;
}

my @sitepp;
for (my $i=0; $i<$totnsite; $i++)	{
	@sitepp[$i] = 0;
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
	$pp[$i] /= $n-$burnin;
	$gene2pp{$genename[$i]} = $pp[$i];
}

open (POSWFILE, $basename.'.sitepp') or die "input error: sitepp file\n";

$n = 0;
foreach my $line (<POSWFILE>)	{
	$n++;
	if ($n > $burnin)	{
        chomp $line;
        my @a = split(/\t+/,$line);
        my $m = (@a);
        if ($m != $totnsite+$Ngene) {
            die "error: non matching number of fields in sitepp file\n";
        }
        my $k = 0;
        my $l = 0;
        for (my $i=0; $i<$Ngene; $i++)	{
            my $name = $a[$k];
            if (! exists $gene2nsite{$name})    {
                die "error: did not find $name in gene list\n";
            }
            $k++;
            if ($l != $gene2offset{$name})  {
                die "non matching offset: $name\t$l\t$gene2offset{$name}\n";
            }
            for (my $j=0; $j<$gene2nsite{$name}; $j++)	{
                $sitepp[$l] += $a[$k];
                $l++;
                $k++;
            }
        }
        if ($l != $totnsite)    {
            die " non matching number of sites\n";
        }
    }
}

for (my $i=0; $i<$totnsite; $i++)	{
	$sitepp[$i] /= $n-$burnin;
}

open (OUTFILE, '>'.$outfile) or die "output error\n";

print OUTFILE "fdr\tgene\tpp\tnsite\tsites\n";
my $totpp = 0;
my $count = 0;
foreach my $gene (sort {$gene2pp{$b} <=> $gene2pp{$a}} keys %gene2pp)	{
	$count++;
	$totpp += $gene2pp{$gene};
	my $fdr = 100 - int(100 * $totpp / $count);
    my $pp = int(100*$gene2pp{$gene});
	print OUTFILE "$fdr\t$gene\t$pp\t";
    my $offset = $gene2offset{$gene};
    my $count = 0;
    for (my $j=0; $j<$gene2nsite{$gene}; $j++)	{
        if ($sitepp[$offset+$j] > $cutoff)  {
            $count++;
        }
    }
    print OUTFILE "$count\t";
    for (my $j=0; $j<$gene2nsite{$gene}; $j++)	{
        if ($sitepp[$offset+$j] > $cutoff)  {
            my $site=$j+1;
            print OUTFILE "$site\t";
        }
    }
    print OUTFILE "\n";
}
	

