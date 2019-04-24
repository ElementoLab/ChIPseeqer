#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndex(0);
my $a_ref = $ta->getArray();


my $r = $h_ref->{$ARGV[1]};

if (!defined($r)) {
  die "Cannot find gene\n";
}

my @b = @$r;
shift @b;

my $or = Sets::order(\@b, 0);

#foreach my $oo (@$or) {
#  print "$oo\t" . $b[$oo] . "\n";
#}




foreach my $r (@$a_ref) {
  my $n = shift @$r; 
  my @a = ();
  for (my $i=0; $i<@$r; $i++) {
    push @a, $r->[ $or->[$i] ];
  }
  print "$n\t" . join("\t", @a) . "\n";
}

