#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;


my $ta = Table->new;

# load table
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;
print join("\t", @$r) . "\n";
my $h_ref = $ta->getIndex(0);


# load and traverse profile
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
$ta->shift;
$ta->sortbycol(1);

foreach my $r (@$a_ref) {
  if (defined($h_ref->{$r->[0]})) {
    print join("\t", @{ $h_ref->{$r->[0]} }) . "\n";
  }
}


