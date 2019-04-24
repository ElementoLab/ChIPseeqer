#!/usr/bin/perl
use lib qw(/home/olly/PERL_MODULES);
use Sets;
use strict;

my $min = (defined($ARGV[1])?$ARGV[1]:0);

open IN, $ARGV[0];
my @COL = ();
my $l = <IN>;
while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    my $n = scalar(@a);
    for (my $i=0; $i<$n; $i++) {
	$COL[$i] ++ if ($a[$i] ne "");
    }
    # print $l if (!Sets::in_array($i, @ARGV));
    # $i ++;
}
close IN;

open IN, $ARGV[0];
while (my $l = <IN>) {    
    chomp $l;
    my @a = split /\t/, $l, -1;
    my @b = ();
    for (my $i=0; $i<scalar(@a); $i++) {
	push @b, $a[$i] if ($COL[$i] > $min);
    }
    print join("\t", @b); print "\n";
}
close IN;



