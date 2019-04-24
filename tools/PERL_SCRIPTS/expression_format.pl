#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;
print     Sets::jointab($r);

foreach my $r (@$a_ref) {
  my $n = shift @$r;
  foreach my $s (@$r) {
    $s = sprintf("%3.2f", $s);
  }
  print "$n\t" . join("\t", @$r) . "\n";_
}

