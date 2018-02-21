use strict;

my $poswfile = shift;
my $Nsite = shift;
my $burnin = shift;
my $outfile = shift;

open (POSWFILE, $poswfile) or die "input error for .posw file: $poswfile\n";
open (OUTFILE, '>'.$outfile) or die "output error: $outfile\n";

my @pp;
for (my $i=0; $i<$Nsite; $i++)	{
	@pp[$i] = 0;
}

my $n = 0;
foreach my $line (<POSWFILE>)	{
	$n++;
	if ($n > $burnin)	{
        chomp $line;
        my @a = split('\t',$line);
        my $nsite = (@a);
        if ($nsite != $Nsite)	{
            die "error when reading .posw file: non matching number of genes: $nsite vs $Nsite\n";
        }
        for (my $i=0; $i<$Nsite; $i++)	{
            if ($a[$i] > 0)	{
                $pp[$i] ++;
            }
        }
	}
}


my %site2pp;
for (my $i=0; $i<$Nsite; $i++)	{
	@pp[$i] /= $n-$burnin;
	$site2pp{$i} = $pp[$i];
}

my $totpp = 0;
my $count = 0;
foreach my $site (sort {$site2pp{$b} <=> $site2pp{$a}} keys %site2pp)	{
	$count++;
	$totpp += $site2pp{$site};
	my $fdr = int(100 * $totpp / $count);
	print OUTFILE "$fdr\t$site\t$site2pp{$site}\n";
}
	


