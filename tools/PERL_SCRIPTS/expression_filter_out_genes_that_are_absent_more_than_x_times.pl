#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;

print join("\t", @$r) . "\n";

foreach my $r (@$a_ref) {
  my $n = shift @$r;
  my $cnt = 0;
  foreach $s (@$r) {
    if ($s eq "") {
      $cnt ++;
    }
  }

  if ($cnt <= $ARGV[1]) {
    print "$n\t" . join("\t", @$r) . "\n";
  }

}

