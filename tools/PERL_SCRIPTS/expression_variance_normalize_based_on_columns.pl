#!/usr/bin/perl

#
# columns number should take into account row names (col 0)
#

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

open IN, $ARGV[0];

use strict;



# index columns to use
my %H = ();
for (my $i=1; $i<@ARGV; $i++) {
  #print STDERR "Norm using $ARGV[$i]\n";
  $H{ $ARGV[$i] } = 1;
}




my %H1 = ();
my %H2 = ();
my $l = <IN>;
print $l;
while (my $l = <IN>) {
  
  chomp $l;
  
  my @a = split /\t/, $l, -1;
  
  # remove ID
  
  
  my @d = ();
  my $i = 0;
  foreach my $r (@a) {

    if (defined($H{$i}) && ($r ne "") && ($r ne 'NaN') && ($r ne 'nan')) {
      push @d, $r;
    }

    $i++;
  }


  # average 
  my $avg = Sets::average(\@d);
  my $std = Sets::stddev (\@d); 

  my @b   = ();
  my $c   = shift @a;

  foreach my $r (@a) {
    if (($r ne "") && ($r ne 'NaN') && ($r ne 'nan')) {

      $r = ( $r - $avg );
      if ($r != 0.0) {
        $r = $r / $std;
      }

      push @b, sprintf("%4.3f", $r);

    } else {

      $r = "";
      push @b, $r;

    }

  }
  
  print $c . "\t" . join("\t", @b) . "\n";
  
}



   
