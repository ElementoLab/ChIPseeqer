#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;

open IN, $ARGV[0];
my $l = <IN>;
#print "RefSeq\tCB/NB log(fold-down)\tLY1/NB log(fold-down)\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  #print "$l\n";  
  
  my $LY1CB = Sets::log2( ($a[3] + 1) / ($a[2] + 1));
  my $LY7CB = Sets::log2( ($a[4] + 1) / ($a[2] + 1));
  my $LY1NB = Sets::log2( ($a[3] + 1) / ($a[1] + 1));
  my $LY7NB = Sets::log2( ($a[4] + 1) / ($a[1] + 1));
  
  if ((abs($LY1CB) > 1) || (abs($LY7CB) > 1)) {
    print "$l\n";
  }
}
close IN;

