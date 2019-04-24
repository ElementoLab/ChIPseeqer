#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Table;
use Sets;

if (@ARGV == 0) {
  die "Args: matrix center[0/1] scale[0/1]\n";
}

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r     = shift @$a_ref;
print     Sets::jointab($r);

my $numlines = @$a_ref;

my $n = @$r;
my @COLS = ();

for (my $i=0; $i<$n; $i++) {
  my $ac = $ta->getColumn($i);
  
  if ($i == 0) {
    push @COLS, $ac;
  } else {
    push @COLS, Sets::center_scale_array($ac, $ARGV[1], $ARGV[2]);
  }
    
}

for (my $i=0; $i<$numlines; $i++) {

  print "$COLS[0][$i]";
  for (my $j=1; $j<$n; $j++) {
    printf("\t%3.2f", $COLS[$j][$i]);
  }
  print "\n";

}



