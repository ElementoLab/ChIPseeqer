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

#print STDERR join("\t", @g1) . "\n";
#print STDERR join("\t", @g2) . "\n";


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

  my $t = Sets::avg_diff(\@g1, \@g2);

  if ($t > Sets::log2(2) && Sets::all_lower(\@g1, \@g2)) {
  #if ($t > Sets::log2(3)) {
    print "$n\t2\n";
  } elsif ($t < Sets::log2(1/2) && Sets::all_greater(\@g1, \@g2)) {
    print "$n\t1\n";
  } else {
    print "$n\t0\n";
  }
  #print "$n\t$t\n";

}

