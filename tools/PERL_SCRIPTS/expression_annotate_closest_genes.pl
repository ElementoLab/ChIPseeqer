#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndex(0);

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray;


foreach my $r (@$a_ref) {
  print "$r->[0]\t$h_ref->{$r->[0]}->[1]\t$h_ref->{$r->[0]}->[2]";

  #print "$r->[0]\t$r->[1]\t$h_ref->{$r->[0]}->[1]";
  my $s = shift @{$h_ref->{$r->[0]}}; shift @{$h_ref->{$r->[0]}};
  if (defined($ARGV[2])) {
    #print "\n    ";
    #print join("\n    ", @{ $h_ref->{$s} } );
  }
  print "\n";
}

