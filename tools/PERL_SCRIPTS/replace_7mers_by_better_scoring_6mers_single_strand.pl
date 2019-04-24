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

my $limit = shift @ARGV;
my $file7 = shift @ARGV;

my $verbose = 0;

#   load the top X 7mers
$t->setLimit($limit);
$t->loadFile($file7);
my $h_ref7 = $t->getIndex(0);
my $a_ref7 = $t->getArray;
$t->setLimit(undef);

#   get the best score
my $lowest_score = $a_ref7->[ scalar(@$a_ref7) - 1 ]->[4];

my $size1        = length($a_ref7->[0]->[0]);


#
#   load the top kmers
#
while (my $filek = shift @ARGV) {
 
    $t->loadFile($filek);
    
    print "Loading $filek\n" if ($verbose == 1);

    my $a_refk_tmp = $t->getArray;
    my $h_refk = ();
    
    #
    #  keep only the top k-mers
    #
    my $a_refk = [];
    foreach my $r (@$a_refk_tmp) {
	last if (sup($lowest_score, $r->[4]));
	$h_refk->{$r->[0]} = $r;
	push @$a_refk, $r;
    }

    my $size2        = length($a_refk->[0]->[0]);

   
    my @OUT = ();

    #
    # SPECIALIZE
    #
    foreach my $r (@$a_ref7) {

	# best candidate is current kmer
	my $bg  = $r;
	
	# initial seq
	my $s = $r->[0];
	
	#
	# look for the best specialization (by adding a nt on both side)
	#
	#foreach my $nt (("A", "C", "T", "G")) {
	    
	my $ss  = substr($s, 1, length($s)-1);

	if ($h_refk->{$ss} && (sup($h_refk->{$ss}->[4], $bg->[4]))) {
	    $bg = $h_refk->{$ss} unless (defined($h_refk->{$ss}->[5]) && ($h_refk->{$ss}->[5] < 10)); 
	}
	
	$ss  = substr($s, 0, length($s)-1); #$nt . $s;
	if ($h_refk->{$ss} && (sup($h_refk->{$ss}->[4], $bg->[4]))) {
	    $bg = $h_refk->{$ss} unless (defined($h_refk->{$ss}->[5]) && ($h_refk->{$ss}->[5] < 10)); 
	} 
	    
	    
	
    
    
	my @aa = ($bg->[0], $bg->[1], $bg->[2], $bg->[3], $bg->[4]);
    
	push @OUT, \@aa;

    }
    
    
    #
    # ADD K-MERS WITHOUT SUBSTRINGS
    #
    #foreach my $r (@$a_refk) {
	

	#next if (defined($r->[5]) && ($r->[5] < 10)); 

	#
	# is there a substring in the refx list
	#
	#my $yes = undef;
	#foreach my $s (@$a_ref7) {
	    # if the 7-mer is a substring of the k-mer OR if the k-mer is a substring of the 7-mer
	    
    #if (Sets::isSubstringOf_1S($s->[0], $r->[0]) || Sets::isSubstringOf_1S($r->[0], $s->[0])) {

				

    #$yes = $s->[0];

		

		#print "include $r->[0]
#}

	#}
	
	
	
	#if (!defined($yes)) {
	    #print "include $r->[0] because they have no substring\n";
	    
	 #   push @OUT, $r;
	#} 
    #}
    
    
    
    #
    #  re-sort the kmers
    #
    @OUT = reverse (sort { main::cmpk($a->[4], $b->[4]) } @OUT);


    #
    #  remove the duplicates
    #
    my @a_new = ();
    my %H = ();
    foreach my $r (@OUT) {
	if (!defined($H{$r->[0]})) {
	    push @a_new, $r;
	} 
	
	$H{$r->[0]} = 1;
    }

    #
    #  remove the lower scoring superstrings
    # 
    
#my %H = ();    
#    my $n = scalar(@a_new);
#    for (my $i=0; $i<$n-1; $i++) {
#	for (my $j=$i+1; $j<$n; $j++) {
#	    if ((length($a_new[$j]->[0]) == length($a_new[$i]->[0])+1) && Sets::isSubstringOf_1S($a_new[$i]->[0], $a_new[$j]->[0])) {
#		$H{$j} = 1;
#	    } 
#	}
#    }
#    my @a_new2 = ();
#    for (my $i=0; $i<$n; $i++) {
#	if (!defined($H{$i})) {
#	    push @a_new2, $a_new[$i];
#	} 
#    }



    #
    #  recreate the list
    #
    
    $a_ref7 = \@a_new;

}


foreach my $r (@$a_ref7) {
    print join("\t", @$r); print "\n";

}

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
