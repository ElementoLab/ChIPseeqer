#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use FileHandle;

use strict;


my @MATRIX = ();
my @COLS   = ();
my @ROWS   = ();

my $cnt_cols = 0;

foreach my $f (@ARGV) {
  
  #print STDERR "Opening $f\n";

  open IN, $f;
  #my $l = <IN>; chomp $l;
  #my @b = split /\t/, $l, -1;

  my $ff = $f;
  $f =~ s/^.+genelist\_//;
  $f =~ s/\_pf\_.+$//;
  push @COLS, $f;

  my $cnt_rows = 0;
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l;
    push @ROWS, $a[0] if ($cnt_cols == 0); 
    $MATRIX[ $cnt_rows ][ $cnt_cols ] = $a[1];
    $cnt_rows ++;
  }
  close IN;

  $cnt_cols ++;
}


foreach my $c (@COLS) {
  print "\t$c";
}
print "\n";

for (my $i=0; $i<@ROWS; $i++) {
  print "$ROWS[$i]";
  for (my $j=0; $j<@COLS; $j++) {
    print "\t" . $MATRIX[$i][$j];
  }
  print "\n";
}



