#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $f = shift @ARGV;

my $ta = Table->new;
$ta->loadFile($f);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r) . "\n";


# go thru @ARGV, assign new index
my %H = ();
my $cnt = 0;
foreach my $r (@ARGV) {
  $H{$r} = $cnt++;
}


foreach my $r (@$a_ref) {
  if (Sets::in_array($r->[1], @ARGV)) {
    print "$r->[0]\t$H{$r->[1]}\n";
  }
}

