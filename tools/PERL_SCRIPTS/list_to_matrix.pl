#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %rownames = ();
my @rowindex = ();
my $cnt_rows = 0;

my %colnames = ();
my @colindex = ();
my $cnt_cols = 0;

my @MATRIX   = ();
foreach my $r (@$a_ref) {

  if (!defined($colnames{$r->[0]})) {
    $colnames{$r->[0]}   = $cnt_cols;
    $colindex[$cnt_cols] = $r->[0];
    $cnt_cols++;
  }

  if (!defined($rownames{$r->[1]})) {
    $rownames{$r->[1]} = $cnt_rows;
    $rowindex[$cnt_rows] = $r->[1];
    $cnt_rows++;
  }

  $MATRIX[ $rownames{$r->[1]} ][ $colnames{$r->[0]} ] = $r->[2];



}

my $m = scalar(keys(%colnames));
my $n = scalar(keys(%rownames));

for (my $i=0; $i<$m; $i++) {
  print "\t$colindex[$i]";
}
print "\n";
for (my $i=0; $i<$n; $i++) {
  print "$rowindex[$i]";
  for (my $j=0; $j<$m; $j++) {
    print "\t$MATRIX[$i][$j]";
  }
  print "\n";

}


