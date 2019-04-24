#
#  take as input a list of genes, and a list of DBs
#
use lib qw(/home/olly/PERL_MODULES);

use Table;
use DataFiles;
use Sequence;
use Getopt::Long;

my $genefile = undef;
my $dbfile   = undef;
my $dir = undef;
GetOptions ('genes=s'          => \$genefile,
	    'dir=s'          => \$dir,
	    'databases=s'      => \$dbfile);



#
#  input list of genes
#

my $a_ref_genes = Sets::readSet($genefile);


#
#  input list of databases
#
my $a_ref_db    = Sets::readSet($dbfile);




my $se = Sequence->new;


if (!$dir) {
    die "Please specify a dir ..\n";
}



my $cnt = 0;
my $a_ref_species = undef;


# 
#  go thru the genes 
#
foreach my $r (@$a_ref_genes) {

    # 
    #  create a file with all sequences
    # 
    open OUT, ">$dir/$r.seq";
    my $m = scalar(@$a_ref_db);
    my $n = 0;
    for (my $i=0; $i<$m; $i++) {

	

	#  set the db
	$se->setBlastDB($a_ref_db->[$i]);

	#  fetch the sequence
	my $seq = $se->getSequenceFromBlastDB($r, 0, 0);
	if ($seq) {

            print "Got $seq from $a_ref_db->[$i]\n";
	    my $ii = $i+1;
	    print OUT ">$r" . "_$ii\n$seq\n\n";
	    $n ++;
	} else {

		print "$r not found ..\n";
        }

        
    }
    close OUT;

    #system("cat $dir/$r->[0].seq");

    #
    # if there are too few sequences, then erase the file
    #
    if ($n < 2) {
	unlink "$dir/$r.seq";
	next;
    }

    #
    #  exec clustalw
    #
    my $s_todo = "clustalw $dir/$r.seq -outorder=input";
    
    print "$s_todo\n";
    
    system($s_todo) == 0 or die "Cannot run Clustalw ..\n";
    

}
