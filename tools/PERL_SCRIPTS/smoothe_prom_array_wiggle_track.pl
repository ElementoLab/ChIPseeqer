#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

open IN, $ARGV[0];

my $prev = 0;
my @set  = ();

while (my $l = <IN>) {
  if (($l =~ /^track/) || ($l =~ /^variable/)) {
    $l =~ s/span\=50/span\=10/;
    $prev = 0;
    print $l;
  } else {

    chomp $l;
    my @a = split /\t/, $l, -1;
    
    # if pos is wthin 500nt of the previous one, add to set
    if ($a[0] <= $prev+1000) {
      push @set, \@a;
      $prev = $a[0]+50;

    } else {
      
      # process current set
      if (@set != 0) {
	
	my $st  = $set[0]->[0];
	my $en  = $set[$#set]->[0]+50;
	my $off = $st;

	print STDERR "$st\t$en\t";
	print STDERR scalar(@set) . "\n";
	
	# fill up array
	my @a_int = ();
	my @a_num = ();
	foreach my $s (@set) {
	  for (my $i=$s->[0]; $i<$s->[0]+50; $i++) {
	    $a_int[$i-$off] += $s->[1];
	    $a_num[$i-$off] ++;
	  }
	}

	# average
	for (my $i=0; $i<$en-$off; $i++) {
	  if ($a_num[$i] > 1) {
	    $a_int[$i] /= $a_num[$i];
	  }
	}

	my @a_newint = ();
	for (my $i=0; $i<$en-$off; $i+=10) {
	  my $sum = 0;
	  my $nit = 0;
	  for (my $j=Sets::max(0,$i-500); $j<Sets::min($i+500, $en-$off); $j++) {
	    my $w = ( 500 - abs($i - $j) ) / 500;
	    if ($a_int[$j] ne "") {
	      $sum += $w * $a_int[$j];
	      $nit += $w;
	    }
	  }
	  my $ii = $i + $off;
	  if ($nit > 0) {
	    $sum /= $nit;
	    print "$ii\t$sum\n";
	  }
	}
	
      }
      
      # create new empty set
      @set = ();
      push @set, \@a;
      $prev = $a[0] + 50;
    }
  } # else
}
close IN;

