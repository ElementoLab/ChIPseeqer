#!/usr/bin/perl

use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;
use Getopt::Long;



my $file1     = $ARGV[0];
my $file2     = $ARGV[1];




open IN1, $file1;
open IN2, $file2;


while (($l1 = <IN1>) && ($l2 = <IN2>)) {
    chomp $l1; chomp $l2;
 
    my @a1 = split /\t/, $l1, -1;
    my @a2 = split /\t/, $l2, -1;
    
    my @l = (@a1, @a2);
    
    print join ("\t", @l) . "\n";
}

close IN1;
close IN2;
