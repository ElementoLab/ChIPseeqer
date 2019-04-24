#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

shift @$a_ref;
my $cnt = 1;
my %C  = ();
foreach my $r (@$a_ref) {
  my $t = shift @$r;
  my $x = join(" ", @$r);

  if (!defined($H{$x})) {
    $H{$x} = $cnt ++;
  }

  $C{$x} ++;

  print "$t\t$H{$x}\n";
  
}


foreach my $k (sort(keys(%H))) {
  print STDERR "$H{$k} => $k \t($C{$k} genes)\n";
}

