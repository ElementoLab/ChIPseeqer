#
#  search for exact conservation of a kmer within a set of aligned sequences
#

use lib qw(/home/olly/PERL_MODULES);




use  ClustalW;
use Sets;


die "Usage : thisprg kmer \"files\"" if (scalar(@ARGV) == 0);


my $a_ref_files = Sets::readSet($ARGV[1]);
#my $a_ref_files = Sets::getFiles($ARGV[1]);

my $cnt_conserved = 0;
my $cnt_there     = 0;
foreach my $f (@$a_ref_files) {

    #print "Opening $f\n";
    
    my $cl = ClustalW->new;
    $cl->setFile($f);

    #print $cl->getHTMLAlignment;

    my $a_ref_aln = $cl->getSeqArray;

    #print "got " . scalar(@$a_ref_aln) . "\n";

    my $b = [];
    my $i = 0;
    my @c = ();
    foreach my $seq (@$a_ref_aln) {
	
	#print "A=$seq\tOO\n";
	
	my @a = ();
	
	my $kmer = $ARGV[0];
	#print "Lookin for $kmer in $seq\n";
	while ($seq =~ /$kmer/ig) {
	    my $p = pos($seq) - length($&);
	    #if ($i == 0) { $cnt_there ++; }
	    $c[$i]++;
	    push @a, $p;
	}

	if ((scalar(@a) > 0) && ($i == 0)) {
	    #print "Found in $f ($seq)\n";
	}

	if ($i == 0) {
	    @$b = @a; 
	} else {
	    #print "OV between \n";
	    #print join("\t", @a); print "\n";
	    #print "AND \n";
	    #print join("\t", @$b); print "\n";
	    
	    $b = Sets::getOverlapSet(\@a, $b);
	}
	
	$i++;
   
	#print "a ";
 	#print join("\t", @a); print "\n";
	
	#print "b ";
	#print join("\t", @$b); print "\n";
    }

    my $cnt = 0;
    foreach my $o (@c) {
	if ($o > 0) {
	    $cnt++;
	}
    }
    $cnt_there ++ if ($cnt == scalar(@$a_ref_aln));
    
    if (scalar(@$b) > 0) {
	#print join("\t", @$b); print "\n";
	$cnt_conserved += scalar(@$b);
    }

    #<STDIN>;
}


print "$ARGV[0]\t$cnt_conserved/$cnt_there\n";
