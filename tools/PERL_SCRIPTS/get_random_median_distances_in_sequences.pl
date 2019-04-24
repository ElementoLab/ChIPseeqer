#
#  input : seqlen array, genes in which to pick the random distances
#




use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;

srand;

#
#  load lengths
#
my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndexKV(0,1);


#
#  load interactions
#
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $a (@$a_ref) {

    #my $todo = "/home/olly/PROGRAMS/MOTIFS_INTERACTIONS/motifs_interaction -kmer1 $a->[1] -kmer2 $a->[2] -fasta1 /home/olly/DATA/WORMS/ENSEMBL_ORTHOLOGS/ELE_D_ACTUAL.seq -fasta2 /home/olly/DATA/WORMS/ENSEMBL_ORTHOLOGS/BRI_D_ACTUAL.seq -nbgenes 11292 -out out.txt.eval -dist dist.txt -twostrand 0 > /dev/null";
    my $todo = "/home/olly/PROGRAMS/MOTIFS_INTERACTIONS/motifs_interaction -kmer1 $a->[1] -kmer2 $a->[2] -fasta1 /home/olly/COMMUNICATIONS/PAPERS/MIRNAS/WEBSITE/DATA/MEL_D_ACTUAL.seq -fasta2 /home/olly/COMMUNICATIONS/PAPERS/MIRNAS/WEBSITE/DATA/PSE_D_ACTUAL.seq -nbgenes 11113 -out out.txt.eval -dist dist.txt -twostrand 0 > /dev/null";

    system($todo);
    
    #
    #  load distances
    #
    $ta->loadFile("dist.txt");
    my $a_ref = $ta->getArray();
    my $a_ref_d = $ta->getColumn(1);
    foreach my $r (@$a_ref_d) {
	$r = abs($r);
    }

    my $n = scalar(@$a_ref_d);
    
    my $d = Sets::median($a_ref_d);
    
    
    my $cnt = 0;
    for (my $i=0; $i<10000; $i++) {
	my @D = ();
	foreach my $r (@$a_ref) {
	    my $l = $h_ref->{ $r->[0] };
	    #my $d = $r->[1];
	    
	    #print "$l\t$d\n";
	    
	    my $x = int(rand($l));
	    
	    my $y = int(rand($l));
	    
	    my $dd =  abs($x - $y); 
	    
	    push @D, $dd;
	}
	my $m =  Sets::median(\@D); #print "\n";
	
	$cnt ++ if ($m < $d);
    }
    
    my $f  = sprintf("%3.2f", $cnt/10000);
    my $fc = sprintf("%3.2f", Sets::min(1.0, ($cnt/10000)*64));
    my $cnn = sprintf("%10s", "$cnt/10000");
    my $sig = ($fc<0.05?"**":"");
	
    print "$a->[1]\t$a->[2]\t$n\t$d\t$cnn\tp(D<d)=$f\tpc=$fc$sig\n"; #<STDIN>;


}
