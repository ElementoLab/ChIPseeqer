#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use FileHandle;



my %MATRIX = ();
my @COLS   = ();

my $cnt_cols = 0;
foreach my $f (@ARGV) {
  open IN, $f;
  my $l = <IN>; chomp $l;
  my @b = split /\t/, $l, -1;
  push @COLS, $b[1];

  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l;
    $MATRIX{ $a[0] }{ $cnt_cols } = $a[1];
  }
  close IN;

  $cnt_cols ++;
}


foreach my $c (@COLS) {
  print "\t$c";
}
print "\n";

foreach my $g (keys(%MATRIX)) {
  print "$g";
  for (my $c=0; $c<@COLS; $c++) {
    print "\t" . $MATRIX{$g}{$c};
  }
  print "\n";
}



