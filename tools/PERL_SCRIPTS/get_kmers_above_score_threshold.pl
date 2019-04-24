#!/usr/bin/perl
use lib qw(/home/olly/PERL_MODULES);
use Sets;


#$a_kmers = Sets::readKmers($ARGV[0]);
open IN, $ARGV[0];

while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l;
    if (($a[4] > $ARGV[1]) || ($a[4] =~ /inf/)) {
	
	print join("\t", @a);
	print "\n";
	
    } else {
      exit;
    }
}




