#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r) . "\n";

my %H = ();



foreach my $r (@$a_ref) {
  push @{ $H{ $r->[0] } }, $r->[1] if (!Sets::in_array($r->[1], @{$H{$r->[0]}}));
}

foreach my $g (keys(%H)) {
  if (scalar(@{$H{$g}}) == 1) {
    print "$g\t$H{$g}->[0]\n";
  } else {
    if ((scalar(@{$H{$g}}) == 2) && Sets::in_array("0", @{$H{$g}})) {

      my $c = undef;
      foreach my $k (@{$H{$g}}) {
	if ($k != 0) {
	  $c = $k;
	  last;
	}	
      }
      die "I don't get it\n" if (!defined($c));

      print "$g\t$c\n";
    }
  }	
}
