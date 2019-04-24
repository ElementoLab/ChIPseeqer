#
#  find group of kmers, instead of single one
#

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";





use Sets;
use GenRegexp;
use DataFiles;
use Table;


my $wmin = 3;

my $w = 500;

my %H = ();
while (my $l = <STDIN>) {
    chomp $l;

    my @a = split /\t/, $l, -1;

    push @{ $H{$a[0]} }, $a[1];
}


foreach my $g (sort(keys(%H))) {

    #
    #  eliminate overlapping instances
    #
    my @aa = @{ $H{$g}};
    my $n = scalar(@aa);
    my @N = ();

    #
    #  sort according to position
    #

    #@aa = sort { $a <=> $b } @aa;

    #print join("\t", @aa) . "\n";
    
    for (my $i=0; $i<$n; $i++) {

	my @win = ();
	
	my $p1 = $aa[$i];
	
	# if start of new cluster overlap with a kmer, go to the next one
	if (($i >= 1) && (($p1 - $aa[$i-1]) < 7)) {
	    next;
	}
	
	#  count the number of non-overlapping occurences
	my $cnt = 1;
	
	# remove all the overlapping instances
	for (my $j=$i+1; $j<$n; $j++) {

	    
	    my $p2 = $aa[$j];
	    my $df  = $p2 - $p1;
	    my $dp  = $p2 - $aa[$j-1];

	    if (($dp >= 3) && ($df <= $w)) {
		$cnt ++;
	    }

	    if ($df > $w) {
		last;
	    }
	    
	}
	
	
	if ($cnt > 0) {
	    
	    $COUNT[ $cnt ] ++;
	    
	    #print "Add cluster of $cnt\n"; 
	}
    }

    #print "Remove cluster of 1 ..\n";
    $COUNT[ 1 ] --;   # remove the terminal cluster ?
    

}

my $n = scalar(@COUNT);
for (my $i=1; $i<$n; $i++) {
    if (!defined($COUNT[$i])) {
	$COUNT[$i] = 0;
    }
    print "$i\t$COUNT[$i]\n";
}



