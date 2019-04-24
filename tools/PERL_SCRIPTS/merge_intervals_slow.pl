#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
my %CHR = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  push @{$CHR{$a[0]}}, \@a;
  
}
close IN;

# chr

foreach my $c (keys(%CHR)) {
 
  my $a_ref_ints = $CHR{$c};
  

  # sort by 
  my @a_int_sorted = sort { $a->[1] <=> $b->[1] } @$a_ref_ints;
  my $numreg = @a_int_sorted;
  
  # 
  my @OUT = ();
  for (my $i=0; $i<$numreg; $i++) {
    $OUT[$i] = 1;
  }

  # compare each reg to all others
  for (my $i=0; $i<$numreg-1; $i++) {
    
    # next if region already deleted
    if ($OUT[$i] == 0) {
      next;
    }
    
    # extend $i as much as possible
    for (my $j=$i+1; $j<$numreg; $j++) {
      # next if region already deleted
      if ($OUT[$j] == 0) {
	next;
      }

      if (Sets::sequencesOverlap($a_int_sorted[$i]->[1], $a_int_sorted[$i]->[2],
				 $a_int_sorted[$j]->[1], $a_int_sorted[$j]->[2])) {
	# extend region $i
	$a_int_sorted[$i]->[1] = Sets::min($a_int_sorted[$i]->[1], $a_int_sorted[$j]->[1]);
	$a_int_sorted[$i]->[2] = Sets::max($a_int_sorted[$i]->[2], $a_int_sorted[$j]->[2]);
	
	if ($a_int_sorted[$i]->[3] !~ /$a_int_sorted[$j]->[3]/) {
	  $a_int_sorted[$i]->[3] = $a_int_sorted[$i]->[3] . "," . $a_int_sorted[$j]->[3];
	}
	if ($a_int_sorted[$i]->[4] !~ /$a_int_sorted[$j]->[4]/) {
	  $a_int_sorted[$i]->[4] = $a_int_sorted[$i]->[4] . "," . $a_int_sorted[$j]->[4];
	}
	$a_int_sorted[$i]->[5] = $a_int_sorted[$i]->[5] . "," . $a_int_sorted[$j]->[5];

	$OUT[$j] = 0;
      } # if overlap	   
      else {
	last; # it's over
      }
      
    } # j loop


    print join("\t", @{$a_int_sorted[$i]}) . "\n";
    
  } # i loop
}


