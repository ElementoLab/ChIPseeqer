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

my $t = $ARGV[1];
if (!defined($t) || ($t eq "")) {
  $t = 1.5;
}

#print "GENE\tEXP\n";
foreach my $r (@$a_ref) {
  my $n = shift @$r;
  #print "$n\t" . join("\t", @$r) . "\n";
  foreach my $s (@$r) {
    if ($s > log($t)) {
      $s = 2;
    } elsif ($s <= -log($t)) {
      $s = 0;
    } else {
      $s = 1;
    }    
  }
  print "$n\t" . join("\t", @$r) . "\n";
}

