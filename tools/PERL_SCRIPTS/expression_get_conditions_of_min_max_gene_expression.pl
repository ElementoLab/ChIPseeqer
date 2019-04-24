#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;


my $ta = Table->new;
$ta->loadFile($ARGV[0]);

# matrix
my $a_ref_g = $ta->getArray();
#shift @$a_ref_g;

# hash table
my $h_ref = $ta->getIndex(0);

# get gene
my @pp = @{$h_ref->{$ARGV[1]}};
my $p  = \@pp;
shift @$p;

# sort exp values
my $o  = Sets::order($p);
#print join("\t", @$o) . "\n"; 

# get n smallest values
my @c1 = ();
my $n  = int(0.5 + $ARGV[2] * @$o);
for (my $i=0; $i<$n; $i++) {
  push @c1, $o->[$i];
}

# get n highest
my @c2 = ();
for (my $i=@$o-$n; $i<@$o; $i++) {
  push @c2, $o->[$i];
}
die "change for loop\n" if (@c2 != $n);


foreach my $r (@$a_ref_g) {
  
  my $n = shift @$r;
  print "$n";
  foreach my $i (@c1) {
    print "\t" . $r->[$i];
  }
  foreach my $i (@c2) {
    print "\t" . $r->[$i];
  }
  print "\n";

}


print STDERR "Got 2 x $n conditions\n";
