#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use FileHandle;



my %MATRIX = ();
my @COLS   = ();

my $cnt_cols = 0;

foreach my $f (@ARGV) {
  
  print STDERR "Opening $f\n";

  open IN, $f;
  #my $l = <IN>; chomp $l;
  #my @b = split /\t/, $l, -1;

  push @COLS, $f;

  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l;
    $MATRIX{ $a[0] }{ $f } = $a[1];
  }
  close IN;

  $cnt_cols ++;
}


foreach my $c (@COLS) {
  print "\t$c";
}
print "\n";

foreach my $g (keys(%MATRIX)) {

  my @v  = ("$g");
  
  my $missing = 0;
  foreach my $c (@COLS) { 
    if (defined($MATRIX{$g}{$c})) {
      push @v, $MATRIX{$g}{$c};
    } else {
      $missing = 1;
    }
  }
  
  if ($missing == 0) {
    print join("\t", @v) . "\n";
  }
}



