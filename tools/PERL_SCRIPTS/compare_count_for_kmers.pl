#!/usr/bin/perl

# prg REAL BKG


# first read bacground
open IN2, $ARGV[1] or die "Cannot open file 2\n";
my %bkg = ();
while (my $l = <IN2>) {
    chomp $l;
    my @a = split /\t/, $l;
    $bkg{ $a[0] } = $a[1];
}
close IN2;


open IN1, $ARGV[0] or die "Cannot open file 1\n";
my @KMERS = ();
while (my $l = <IN1>) {
    chomp $l;
    my @a = split /\t/, $l;
    print "$a[0]\t$a[1]\t" . $bkg{$a[0]} . "\t" . sprintf("%3.2f", $a[1]/$bkg{$a[0]}) . "\n"; #$bkg{ $a[0] }->[2] = $a[1]; 
}
close IN1;




