#
#  find group of kmers, instead of single one
#

#  input : entity, chromosome, position

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Sets;
use Table;
use strict;

if (!$ARGV[0] && !$ARGV[1]) {
  die "please input a position file and a window size\n";

}


my $ov = 0;
my $wi = $ARGV[1];
my $nm = 2;


#my $a_ref_pos_1 = Sets::readSet($ARGV[0]);

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref_pos_1 = $ta->getArray();


#
#  get sublists for each chromosome
#
my %CHR = ();
foreach my $p (@$a_ref_pos_1) {

    push @{ $CHR{ $p->[1] } }, $p;

}


foreach my $c (keys(%CHR)) {
    
    

    #
    #  sort all the positions of a single chromosome
    #
    my @a         = @{ $CHR{ $c } };
    my @a         = sort { $a->[2] <=> $b->[2] } @a;

    my $a_ref_pos = \@a;
        
    
    my @P = ();
    if ($ov > 0) {
	
	my $n = scalar(@$a_ref_pos);
	my @N = ();
	#	
	# detect all the overlapping instances
	#
	for (my $i=0; $i<$n-1; $i++) {
	    for (my $j=$i+1; $j<$n; $j++) {
		my $p1 = $a_ref_pos->[$i]; 
		my $p2 = $a_ref_pos->[$j];
		if (abs($p1-$p2) < $ov) {
		    $N[$j] = 1;
		}
	    }
	}
		
	#
	#  get only non-overlapping positions
	#
	@P = ();
	for (my $i=0; $i<$n; $i++) {
	    if ($N[$i] != 1) {
		push @P, $a_ref_pos->[$i];
	    }
	}
	
    } else {
	@P = @$a_ref_pos;
    }
    
    
    
    
    #
    #  calculate all clusters (of size $wi nt)
    #
    my $n = scalar(@P);
    my @ENHANCERS = ();
    for (my $i=0; $i<$n-1; $i++) {
	my @enhancer = ();
	push @enhancer, $P[$i];
	for (my $j=$i+1; $j<$n; $j++) {
	    if (abs($P[$i]->[2] - $P[$j]->[2]) <= $wi) {
		push @enhancer, $P[$j];
	    } else { 
		last;
	    }
	    
	}
	
	
	if (scalar(@enhancer) >= $nm) {
	    #push @ENHANCERS, \@enhancer;
	    print "Cluster\n";
	    foreach my $e (@enhancer) {
		print join("\t", @$e); print "\n";
	    }
	}
    }
	
}
