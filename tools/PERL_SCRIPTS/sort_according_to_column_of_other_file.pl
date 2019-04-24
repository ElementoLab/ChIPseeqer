#!/usr/bin/perl

#
#  usage perl ~/PERL_MODULES/SCRIPTS/sort_according_to_column_of_other_file.pl X file1 col_to_sort file_to_sort
#

use lib qw(/home/elemento/PERL_MODULES);
use Table;
use Sets;

Sets::cmdLine(1, "X file1 col_to_sort file_to_sort");

my $cc1 = $ARGV[0];
my @b = ();
open IN, $ARGV[1];
while ($l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    push @b, $a[$cc1];
}
close IN;



my $cc2 = $ARGV[2];

my $ta = Table->new;
$ta->loadFile($ARGV[3]);
my $h_ref = $ta->getIndex($cc2);

foreach my $r (@b) {
    print join("\t", @{ $h_ref->{$r} }); print "\n";
}




