#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;
print "$r->[0]\tt-value\n";
my $n1 = $ARGV[1];
my $n2 = $ARGV[2];


my $n = shift @$r;

my @g1 = ();
for (my $i=0; $i<$n1; $i++) {    
  push @g1, $r->[$i];
}

my @g2 = ();
for (my $i=$n1; $i<$n1+$n2; $i++) {
  push @g2, $r->[$i];
}

print STDERR join("\t", @g1) . "\n";
print STDERR join("\t", @g2) . "\n";


foreach my $r (@$a_ref) {
  my $n = shift @$r;

  my @g1 = ();
  for (my $i=0; $i<$n1; $i++) {    
    push @g1, $r->[$i];
  }

  my @g2 = ();
  for (my $i=$n1; $i<$n1+$n2; $i++) {
    if ($r->[$i] eq "") {
      print STDERR "Careful, undef value at col $i\n";
	die "Pb\n";
    }
    push @g2, $r->[$i];
  }

  my $t  = Sets::t_statistic(\@g1, \@g2);

  my $a1 = Sets::average(\@g1);
  my $a2 = Sets::average(\@g2);

  print "$n\t$t\t$a1\t$a2\n";

}

