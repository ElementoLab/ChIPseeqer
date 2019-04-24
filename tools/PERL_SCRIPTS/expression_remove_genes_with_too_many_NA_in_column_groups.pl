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
print join("\t", @$r) . "\n";

my $n1 = $ARGV[1];
my $n2 = $ARGV[2];

my $t  = 3;
if ($ARGV[3] ne "") {
  $t = $ARGV[3];
}

foreach my $r (@$a_ref) {
  my $n = shift @$r;
  
  my $cntNA1 = 0;
  for (my $i=0; $i<$n1; $i++) {    
    $cntNA1++ if ($r->[$i] eq "NA");
  }
  
  my $cntNA2 = 0;
  for (my $i=$n1; $i<$n1+$n2; $i++) {
    $cntNA2++ if ($r->[$i] eq "NA");
  }
  
  next if ( $n1 - $cntNA1 < $t);
  next if ( $n2 - $cntNA2 < $t);

  print "$n\t" . join("\t", @$r) . "\n"; 

}

