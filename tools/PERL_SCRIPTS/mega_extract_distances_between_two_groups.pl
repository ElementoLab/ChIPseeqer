#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use strict;

open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";

my @hyperm = split /\,/, $ARGV[3];

my $l = undef;
my @dist   = ();
my %codes  = ();

my $in = 0;
while ($l = <IN>) { 
  chomp $l;
  $l =~ s/\r//g;
  if ($l =~ /^\[/) {
    $in = 1;
  }

  if ($in == 1) {
    if ($l =~ /^\[\ *(\d+)\]\ \#(.+?)$/) {

      #print "$2 => $1\n\n";
      $codes{$2} = $1;

    } elsif ($l =~ /^\[\ *(\d+)\]\ \ (0.+?)$/) {

      my $lab  = $1;
      my $rest = $2;

      my @a = split /\ +/, $rest;
      #shift @a;

      for (my $i=0; $i<@a; $i++) {
	$dist[$lab][$i+1] = $a[$i];
	$dist[$i+1][$lab] = $a[$i];
      }

      #print "$lab\t" . scalar(@a) . "\t$l\n";

    }
  }

}


close IN;

my @seqs = keys(%codes);

for (my $i=1; $i<@seqs; $i++) { 
  my $s1 = $seqs[$i];

  next if ($s1 !~ /$ARGV[1]/);
  next if (Sets::in_array($s1, @hyperm));

  for (my $j=1; $j<@seqs; $j++) { 
    my $s2 = $seqs[$j];
      
    next if ($s2 !~ /$ARGV[2]/);
    next if (Sets::in_array($s2, @hyperm));

    my $d = $dist[ $codes{$s1} ][ $codes{$s2} ];
    

    print "$s1\t$s2\t$d\n";
    
    if (!defined($d)) {
      die "Problem with $s1 and $s2\n";
    }

  }

}
