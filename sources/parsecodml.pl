
use strict;

my $infile = shift;
my $outfile = shift;

open(INFILE, $infile) or die "input error\n";
open(OUTFILE, '>'.$outfile) or die "output error\n";

my %gene2delta;

while (my $line = <INFILE>) {

    chomp $line;
    if ($line =~ /^reduced(.+)\.codeml.+\s+(\-\d+\.\d+)\s/g)    {
        my $name = $1;
        my $lnl1 = $2;
        $line = <INFILE>;
        chomp $line;
        if ($line =~ /^reduced(.+)\.codeml.+\s+(\-\d+\.\d+)\s/g)    {
            if ($name ne $1)    {
                die "error: not twice the same gene name\n";
            }
            my $lnl2 = $2;
            my $delta = $lnl2 - $lnl1;
            $gene2delta{$name} = $delta;
        }
        else    {
            die "error when parsing line2: $line\n";
        }
    }
    else    {
        die "error when parsing line1: $line\n";
    }
}

foreach my $gene (sort {$gene2delta{$b} <=> $gene2delta{$a}} keys %gene2delta)	{
    my $delta = $gene2delta{$gene};
    print OUTFILE "$gene\t$delta\n";
}
