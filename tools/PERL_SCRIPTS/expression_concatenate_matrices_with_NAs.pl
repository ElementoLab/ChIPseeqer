#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use FileHandle;
use strict;


my %MATRIX   = ();
my @COLS     = ();
my $cntcols = 0;
my @ROWS     = ();
my %HROWS    = ();

foreach my $f (@ARGV) {
  
  print STDERR "STDERR (not an error): Opening $f\n";

  open IN, $f;

  # first row is special
  my $l = <IN>; chomp $l;
  my @b = split /\t/, $l, -1;

  my $n = shift @b;  # ignore
  my $m = @b;
  push @COLS, @b;



  while (my $l = <IN>) {

    chomp $l;
    my @a = split /\t/, $l, -1;
    my $n = shift @a;

    for (my $i=0; $i<$m; $i++) {
      $MATRIX{$n}[$cntcols+$i] = $a[$i];
    }
    
  }
  
  $cntcols += $m;
  
  close IN;

}


print "ID";
foreach my $c (@COLS) {
  print "\t$c";
}
print "\n";

my @ROWS = keys(%MATRIX);

foreach my $g (@ROWS) {
  print "$g";
  for (my $i=0; $i<$cntcols; $i++) {
    print "\t" . $MATRIX{$g}->[$i]; 
  }
  print "\n";
}



