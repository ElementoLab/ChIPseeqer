#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use strict;

my $ap2list = "$ENV{HOME}/PROGRAMS/PLASMODIUM/PBM/PWMs/RAW/PBMCORES/AP2_PBM_list";

my @dirs    = @ARGV;

my $a_ref = Sets::readSet($ap2list);

foreach my $d (@dirs) {
  print "\t$d";
}
print "\n";

foreach my $ap2 (@$a_ref) {
  print "$ap2";
  foreach my $d (@dirs) {
    
    my %GENES = ();
    open IN, "$d/$ap2.matches.txt";
    while (my $l = <IN>) {
      chomp $l;
      my @a = split /\t/, $l, -1;
      $GENES{$a[0]}++;
    }
    close IN;
    
    my $cnt = 0;
    foreach my $g (keys(%GENES)) {
      $cnt += $GENES{$g};
    }
    $cnt /= scalar(keys(%GENES));

    print "\t" . scalar(keys(%GENES)) . sprintf("(%3.1f)", $cnt);
  }

  print "\n";
  
  
  
}
