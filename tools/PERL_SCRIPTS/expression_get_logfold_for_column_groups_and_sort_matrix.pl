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
#print "$r->[0]\tt-value\n";
print join("\t", @$r) . "\n";
my $n1 = $ARGV[1];
my $n2 = $ARGV[2];


my @ts = ();
my $i  = 0;
foreach my $r (@$a_ref) {
  #my $n = shift @$r;

  my @g1 = ();
  for (my $i=0; $i<$n1; $i++) {
    push @g1, $r->[$i+1];
  }

  my @g2 = ();
  for (my $i=$n1; $i<$n1+$n2; $i++) {
    if ($r->[$i+1] eq "") {
      print STDERR "Careful, undef value at col $i\n";
    }
    push @g2, $r->[$i+1];
  }

  #my $t = Sets::t_statistic(\@g1, \@g2);
  my $t = Sets::log2( (Sets::average(\@g1) + 1) / ( Sets::average(\@g2) + 1) );
  my @a_t = ($i, $t);

  push @ts, \@a_t;
  
  #print "$n\t$t\n";

  $i ++;
}

@ts = sort { $b->[1] <=> $a->[1] } @ts;

foreach my $r (@ts) {
  print join("\t", @{ $a_ref->[$r->[0]] }) . "\n";
}
