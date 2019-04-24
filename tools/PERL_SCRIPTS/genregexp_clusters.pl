#
#  find group of kmers, instead of single one
#

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";




use ClustalW;
use Sets;
use GenRegexp;
use DataFiles;
use Table;
use Drosophila;

if (!$ARGV[0] && !$ARGV[1]) {
  die "please input a fasta file and a list of motifs\n";

}

my $dr = Drosophila->new;
my $ov = 4;
my $fasta = $ARGV[0];
my $a_ref_block = Sets::readSet($ARGV[1]);

my %H = ();
my $cnt = 1;
foreach my $k (@$a_ref_block) {
    
    my $ge = GenRegexp->new;
    
    if ($ARGV[2] eq '0') {
	$ge->setTwoStrands(0);
    } 
    
    $ge->calc($k, $fasta);
    
    my $ta = $ge->getTable; #getNbMatches(); 
    
    my $a_ref = $ta->getArray;

    foreach my $r (@$a_ref) {

	#print "$r->[0]\n";
	push @{$H{$r->[0]} }, "$cnt($r->[1], $r->[2])";
    }
    
    #print "$k\t$c\n";

    $cnt ++;
}


my @P = ();
foreach my $g (sort(keys(%H))) {

    #
    #  eliminate overlapping instances
    #
    my @a = @{ $H{$g}};
    my $n = scalar(@a);
    my @N = ();

    for (my $i=0; $i<$n-1; $i++) {
	
	# remove all the overlapping instances
	for (my $j=$i+1; $j<$n; $j++) {

	    my ($p1) = $a[$i] =~ /\d+\((\d+)/;
	    my ($p2) = $a[$j] =~ /\d+\((\d+)/;
	    

	    #print " abs($a[$j]-$a[$i]) < 7 ?\n";
	    
	    
	    
	    if (abs($p1-$p2) < $ov) {
		$N[$j] = 1;
	    }
	    
	}
	
    }
    
    
    my @D = ();
    my $name = $dr->getGeneName($g);
    print "$g\t$name\t";
    my $nbn1 = 0;
    for (my $i=0; $i<$n; $i++) {


	if ($N[$i] != 1) {
	    my ($n1, $p1) = $H{$g}->[$i] =~ /(\d+)\((\d+)/;
	    if ($n1 == 1) {
		$nbn1 ++;
	    }
	    #print $H{$g}->[$i]; print "\t";
	    push @D, $H{$g}->[$i];
	    push @P, $p1;
	}
    }
    
    # what is the maximum number that you 
    #   can put in a 400bp window
    my @pos = ();
    foreach my $o (@D) {
	my @a_tmp = $o =~ /(\d+)\((\d+)\)/;
	push @pos, \@a_tmp;
    }

    @pos = sort { $a->[1] <=> $b->[1] } @pos;
    
    # order them by position
    my $n = scalar(@pos);
    my @ENHANCERS = ();
    for (my $i=0; $i<$n-1; $i++) {
	# start a new enhancer
	my @enhancer = ();
	push @enhancer, $pos[$i];
	for (my $j=$i+1; $j<$n; $j++) {
	    
	    if (abs($pos[$i]->[1] - $pos[$j]->[1]) <= 400) {
		push @enhancer, $pos[$j];
	    }
	    
	}

	if (scalar(@enhancer) >= 3) {
	    push @ENHANCERS, \@enhancer;
	}
	
    }

    my $ne = scalar(@ENHANCERS);
    
    print scalar(@D); print "\t";
    #print "$nbn1\t";
    #print "$ne\t";
    print join("\t", @D); print "\n"; 
    #print "\n";

    
}


Sets::printSet(\@P);
