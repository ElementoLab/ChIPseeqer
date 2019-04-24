use strict;
open IN, $ARGV[0];

my %GENES = ();

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;


  push @{ $GENES{$a[0]}{$a[1]} }, \@a;


}
close IN;

foreach my $g (keys(%GENES)) {
  
  foreach my $m (keys(%{$GENES{$g}})) {
    
    #
    if ($m eq $ARGV[1]) {
      
      my @l = sort { $a->[2] <=> $b->[2] } @{$GENES{$g}{$m}};
      
      my @TD = ($l[0]);

      for (my $i=1; $i<@l; $i++) {
	
	if (@TD == 0) {
	  push @TD, $l[$i];

	} else {
	  
	  my $d = $l[$i]->[2] - $TD[0]->[2];
	  if (($d > 0) && ($d <= 15)) {
	    push @TD, $l[$i];
	    # avg pos
	    my $pos = 0;
	    my $sco = 0;
	    my @occ = ();
	    my @apos = ();
	    foreach my $t (@TD) {
	      $pos += $t->[2];
	      $sco += $t->[4];
	      push @occ, $t->[5];
	      push @apos, $t->[2];
	    }
	    $pos /= scalar(@TD);
	    $sco /= scalar(@TD);
	    $pos = int(0.5+$pos);
	    print "$g\t$m-SPE2\t$pos\t1\t$sco\t" . join("/", @occ) . "\n";
	    #print "$g\t$m-Tandem\t". join("/", @apos) . "//$pos\t1\t$sco\t" . join("/", @occ) . "\n";
	    @TD = (); # empty store
	  } else {
	    
	    # too distant, print store, empty it, add new
	    print join("\t", @{$TD[0]}) . "\n";
	    @TD = ($l[$i]);
	    
	  }
	} # else 
	  
      } # end for
      if (@TD == 1) {
	print join("\t", @{$TD[0]}) . "\n";
      }

    } # if right motif
    
  }
}
