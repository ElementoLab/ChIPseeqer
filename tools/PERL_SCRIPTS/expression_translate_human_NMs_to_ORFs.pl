#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile("/Users/olivier/PROGRAMS/PAGE/PAGE_DATA/ANNOTATIONS/human_go_orf_old/human_go_orf_old_genedesc.txt");
my $h_ref = $ta->getIndex(1);

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray;
my $r = shift @$a_ref;
print join("\t", @$r) . "\n";

foreach my $r (@$a_ref) {
  print $h_ref->{$r->[0]}->[0] . "\t$r->[1]\n"; 

}

