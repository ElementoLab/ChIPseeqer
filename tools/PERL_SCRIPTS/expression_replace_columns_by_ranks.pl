#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $nc = $ta->getNbColumns();

#print "$nc columns\n";

my $a_ref = $ta->getArray();
my $n     = @$a_ref;

for (my $i=1; $i<$nc; $i++) {

  my $a_ref_col   = $ta->getColumn($i);
  shift @$a_ref_col;
  my $a_ref_ranks = Sets::ranks($a_ref_col);

  #for (my $j=0; $j<@$a_ref_col; $j++) {
  #  print "$a_ref_col->[$j]\t$a_ref_ranks->[$j]\n";
  #}
  #<STDIN>;

  #print join("\t", @$a_ref_ranks) . "\n";

  my $j = 0;
  my $l = 0;
  foreach my $r (@$a_ref) {
    if ($l > 0) {
      $r->[$i] = $a_ref_ranks->[$j] / $n;
      $j++;
    }
    $l++;
  }

}


#my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  print join("\t", @$r) . "\n";
}

