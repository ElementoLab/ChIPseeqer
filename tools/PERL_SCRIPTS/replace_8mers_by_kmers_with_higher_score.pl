#
#   starting from the X top 7-mers, try to improve the score by looking at generalized (6mers) or specialized (8mers) kmers
#

use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;
use strict;

my $t = Table->new;


my $size1 = undef;
my $size2 = undef;


#   load the X 8mers


$t->setLimit($ARGV[0]);
$t->loadFile($ARGV[1]);
my $h_ref7 = $t->getIndex(0);
my $a_ref7 = $t->getArray;
$t->setLimit(undef);

#   get the best score
my $lowest_score = $a_ref7->[ scalar(@$a_ref7) - 1 ]->[4];


#   load the top X mers
$t->loadFile($ARGV[2]);
my $a_ref8_tmp = $t->getArray;
my $h_ref8 = ();
my $a_ref8 = [];
foreach my $r (@$a_ref8_tmp) {
    last if (sup($lowest_score, $r->[4]));
    $h_ref8->{$r->[0]} = $r;

    push @$a_ref8, $r;
}


#   get the size of the list 
my $size2        = length($a_ref8->[0]->[0]);

   
my @OUT = ();
# traverse the array, trying to find new members
foreach my $r (@$a_ref7) {
    
    #   get the size 
    my $size1        = length($r->[0]);

    #print "$r->[0] .$r->[4].\n";

    # best candidate
    my $bg  = $r;
	    
    # initial seq
    my $s = $r->[0];
    
    if ($size1 < $size2) {
	
	#print "specialize ..\n";

	# find best specialization (by adding a nt on both side)
	foreach my $nt (("A", "C", "T", "G")) {

	    
	    my $ss  = $s . $nt;
	    my $ssc = Sets::getComplement($ss);
	    if ($h_ref8->{$ss} && (sup($h_ref8->{$ss}->[4], $bg->[4]))) {
		$bg = $h_ref8->{$ss};
	    } 
	    if ($h_ref8->{$ssc} && (sup($h_ref8->{$ssc}->[4], $bg->[4]))) {
		$bg = $h_ref8->{$ssc};
	    } 
	    

	    $ss  = $nt . $s;
	    $ssc = Sets::getComplement($ss);
	    if ($h_ref8->{$ss} && (sup($h_ref8->{$ss}->[4], $bg->[4]))) {
		$bg = $h_ref8->{$ss};
	    } 
	    if ($h_ref8->{$ssc} && (sup($h_ref8->{$ssc}->[4], $bg->[4]))) {
		$bg = $h_ref8->{$ssc};
	    } 


	}

    } elsif ($size1 > $size2) {
	
	#
	#   code to specialize motifs 
	#
	
	my $g   = 1;
	
	#
	#   calc the substring 1
	#
	my $ss  = substr($s, 0, length($s) - $g);
	my $ssc = Sets::getComplement($ss);

	#print "$ss\t$ssc\t$g\n";

	if ($h_ref8->{$ss} && (sup($h_ref8->{$ss}->[4], $bg->[4]))) {
	    $bg = $h_ref8->{$ss};
	} 
	if ($h_ref8->{$ssc} && (sup($h_ref8->{$ssc}->[4], $bg->[4]))) {
	    $bg = $h_ref8->{$ssc};
	} 
		
	#
	#   calc the substring 2
	#
	my $ss  = substr($s, $g, length($s) - $g);
	my $ssc = Sets::getComplement($ss);
	if ($h_ref8->{$ss} && (sup($h_ref8->{$ss}->[4], $bg->[4]))) {
	    $bg = $h_ref8->{$ss};
	} 
	if ($h_ref8->{$ssc} && (sup($h_ref8->{$ssc}->[4], $bg->[4]))) {
	    $bg = $h_ref8->{$ssc};
	} 
	
    }
    


    
    
    #print "$bg->[0]\t$bg->[1]\t$bg->[2]\t$bg->[3]\t$bg->[4]\n";
    
    my @aa = ($bg->[0], $bg->[1], $bg->[2], $bg->[3], $bg->[4]);
    
    push @OUT, \@aa;

}


if ($ARGV[3]) {
    foreach my $r (@$a_ref8) {
	
	# is there a substring in the refx list
	my $yes = undef;
	foreach my $s (@$a_ref7) {
	    
	    #if (length($s->[0]) != length($r->[0])+1) {
	    #    $yes = "NNN";
	    #    next;
	    #}
	    
	    #print "is $s->[0] a substring of $r->[0] ?\n";
	    if (Sets::isSubstringOf_2S($s->[0], $r->[0]) || Sets::isSubstringOf_2S($r->[0], $s->[0])) {
		$yes = $s->[0];
	    }
	}
	
	
	
	if (!defined($yes)) {
	    #print "NO SUBS .." . join("\t", @$r) . "\n";
	    push @OUT, $r;
	} else {
	    #print "$r->[0] has a substring $yes ..\n";
	}
    }
}


@OUT = reverse (sort { main::cmpk($a->[4], $b->[4]) } @OUT);

#
#  remove lower scoring substrings
#
my $n = scalar(@OUT);
for (my $i=0; $i<$n-1; $i++) {
   for (my $j=$i+1; $j<$n; $j++) {
       
       if (Sets::isSubstringOf_2S($OUT[$j]->[0], $OUT[$i]->[0], 1)) {
	   $OUT[$j]->[5] = 'N';
       }
   } 
}

#
# eliminate doubles and display
#
my %H = ();
foreach my $r (@OUT) {
    if (!defined($H{$r->[0]}) && ($r->[5] ne 'N')) {
	print join("\t", @$r) . "\n";
    } 

    $H{$r->[0]} = 1;
}

exit;














#my @OUTNSS = ();







sub cmpk {
    my ($f1, $f2) = @_;

    if ($f1 =~ /inf/) {
	return 1;
    }

    if ($f2 =~ /inf/) {
	return -1;
    }
	
    return ($f1 <=> $f2);
    
}


sub sup {
    my ($f1, $f2) = @_;

    if ($f1 =~ /inf/) {
	return 1;
    }

    if ($f2 =~ /inf/) {
	return 0;
    }
	
    return ($f1 > $f2);
    
}
