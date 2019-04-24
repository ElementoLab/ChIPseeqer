use strict;

my @ls = <STDIN>; chomp @ls;


my $n = scalar(@ls);

for (my $j=0; $j<$n; $j++) {
    print "\t" . substr($ls[$j], 0, 4);
}
print "\n";

for (my $i=0; $i<$n; $i++) {
    print substr($ls[$i], 0, 4);

    for (my $j=0; $j<$n; $j++) {
 
	

	my @seq1 = split //, $ls[$i];
	my @seq2 = split //, $ls[$j];
	
	my $i = 4;
	my $lgtrue = 0;
	my $p = 0;
	while( defined($seq1[$i] ) && defined($seq2[$i]) ) {
	    #print "$seq1[$i] => $seq2[$i]\n";
	    if(($seq1[$i] ne '-') && ($seq2[$i] ne '-')) { 
		$lgtrue++;
		if ($seq1[$i] ne $seq2[$i]) { 
		    $p++; #print " $p";
		}
	    }
	    $i++;
	}
	
	$p = $p / $lgtrue;
	
	if (1-20.*$p/19.<=0) { 
	    print "$p\t-1\n";
	}
	
	#print $p;
	my $k= log (1 - 20*$p/19)*(-19/20);
	print "\t" . sprintf("%4.3f", $k);
    }
    print "\n";
}
