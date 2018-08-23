
use strict;

my $naa = 20;

my $basename = shift;
my $cond = shift;
my $burnin = shift;
my $cutoff = shift;
my $offset = shift;

open (TOGGLEFILE, $basename.'_'.$cond.'.shifttoggle') or die "input error: shifttoggle file\n";

my %siteaapp;
my %sitepp;
my $nsite = 0;
my $n = 0;

foreach my $line (<TOGGLEFILE>)	{
	$n++;
	if ($n > $burnin)	{
        chomp $line;
        my @a = split(/\t+/,$line);
        my $m = (@a);

        if (! $nsite)   {
            $nsite = $m;
            for (my $i=0; $i<$nsite*$naa; $i++)	{
                $siteaapp{$i} = 0;
            }
            for (my $i=0; $i<$nsite; $i++)	{
                $sitepp{$i} = 0;
            }
        }
        else    {
            if ($nsite != $m)   {
                die "non matching number of entries\n";
            }
            for (my $i=0; $i<$nsite*$naa; $i++)	{
                $siteaapp{$i} += $a[$i];
            }
            my $k = 0;
            for (my $i=0; $i<$nsite; $i++)	{
                my $l = 0;
                for (my $m=0; $m<$naa; $m++)    {
                    if ($a[$k]) {
                        $l = 1;
                    }
                    $k++;
                }
                if ($l) {
                    $sitepp{$i}++;
                }
            }
        }
    }
}

for (my $i=0; $i<$nsite; $i++)	{
	$sitepp{$i} /= $n-$burnin;
}

for (my $i=0; $i<$nsite*$naa; $i++)	{
	$siteaapp{$i} /= $n-$burnin;
}

open (OUTFILE, '>'.$basename.'.sortedsitepp');
print OUTFILE "pos\tpp\tfdr";
print OUTFILE "\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY";
print OUTFILE "\n";

my $totpp = 0;
my $ndisc = 0;
foreach my $i (sort {$sitepp{$b} <=> $sitepp{$a}} keys %sitepp) {

    if ($sitepp{$i} > $cutoff)  {
        $ndisc++;
        $totpp += $sitepp{$i};
        my $f = $totpp / $ndisc;
        my $fdr = int(100 * $f);
        my $pp = int(100 * $sitepp{$i});
        my $pos = $i + $offset;
        print OUTFILE "$pos\t$pp\t$fdr";
        for (my $a=0; $a<$naa; $a++)    {
            my $aapp = int(100*$siteaapp{$i*$naa+$a});
            print OUTFILE "\t$aapp";
        }
        print OUTFILE "\n";
    }
}



