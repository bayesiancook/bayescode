use strict;

my $infile = shift;
my $outfile = shift;
my $burnin = shift;

open (INFILE, $infile) or die "input error\n";
open (OUTFILE, '>'.$outfile) or die "output error\n";

my $line = <INFILE>;
chomp $line;
my @a = split('\t', $line);
my $nparam = (@a);

my @mean;
for (my $i=0; $i<$nparam; $i++) {
    $mean[$i] = 0;
}

my $count = 1;
my $tot = 0;
foreach my $line (<INFILE>) {
    chomp $line;
    $count++;
    if ($count > $burnin)   {
        my @a = split('\t',$line);
        for (my $i=0; $i<$nparam; $i++) {
            $mean[$i] += $a[$i];
        }
        $tot++;
    }
}

for (my $i=0; $i<$nparam; $i++) {
    $mean[$i] /= $tot;
    print OUTFILE "$mean[$i]\t";
}
print OUTFILE "\n";

