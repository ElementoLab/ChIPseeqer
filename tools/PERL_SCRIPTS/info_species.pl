#!/usr/bin/perl

use lib qw(/home/olly/PERL_MODULES);

use DataFiles;


my $df = DataFiles->new;

my $r = $df->get($ARGV[0]);

foreach my $k (sort(keys(%$r))) {
    print "$k\t$r->{$k}\n";
}
