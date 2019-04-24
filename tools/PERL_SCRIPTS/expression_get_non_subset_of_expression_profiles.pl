#!/usr/bin/perl

use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use Table;
use strict;

#my $a_ref = Sets::readSet($ARGV[0]);

my $h_ref = Sets::getIndexFromTableColumn($ARGV[0], 0);

#print scalar(keys(%$h_ref));


open IN, $ARGV[1];
if (!defined($ARGV[2])) {
  my $l = <IN>;
  print "$l";
}
while (my $l = <IN>) {
    chomp $l;

    my @a = split /\t/, $l, -1;
    
    if (!defined($h_ref->{ $a[0] })) {
	print "$l\n";
    }

}
close IN;
