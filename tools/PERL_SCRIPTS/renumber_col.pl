#!/usr/bin/perl

use lib qw(/home/elemento/PERL_MODULES);

use Table;
use Sets;





my $cnt = 1;

if (!defined($ARGV[1])) {
  $l1 = <STDIN>;
  print $l1;
}
while ($l1 = <STDIN>) {
    chomp $l1; 
 
    my @a1 = split /\t/, $l1, -1;

    $a1[ $ARGV[0] ] = $cnt;

    print join("\t", @a1); print "\n";

    $cnt++;

}

