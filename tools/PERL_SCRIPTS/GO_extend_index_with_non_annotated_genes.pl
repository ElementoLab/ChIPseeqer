#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getColumn(0);

$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndex(0);


foreach my $r (@$a_ref) {
  if (defined($h_ref->{$r})) {
    print join("\t", @{$h_ref->{$r}}) . "\n";
  } else {
    print "$r\n"
  }
}

