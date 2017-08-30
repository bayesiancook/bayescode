use strict;

my $infile = shift;
my $outfile = shift;

open(INFILE, $infile) or die "input error\n";
open(OUTFILE, '>'.$outfile) or die "output error\n";

my $header = <INFILE>;
chomp $header;
my @a = split('\t',$header);
my $n = (@a);
if ($n != 2)    {
    die "header does not contain 2 fields\n";
}
my $ntaxa = $a[0];
my $nsite = $a[1];

my %tax2seq;
my $ntax = 0;
foreach my $line (<INFILE>) {

    chomp $line;
    my @a = split('\t',$line);
    my $n = (@a);
    if ($n != 2)    {
        die "line does not contain 2 fields\n";
    }
    my $tax = $a[0];
    my $seq = $a[1];
    if ($seq =~ /^[?-X]+$/g)    {
        print "$tax all missing\n";
    }
    else    {
        $ntax ++;
        $tax2seq{$tax} = $seq;
    }
}

print OUTFILE "$ntax\t$nsite\n";
foreach my $tax (keys %tax2seq) {
    print OUTFILE "$tax  $tax2seq{$tax}\n";
}

