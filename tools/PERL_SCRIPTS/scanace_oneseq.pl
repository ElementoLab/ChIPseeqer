use lib qw(/home/olly/PERL_MODULES);
require "ScanACE.pm";
require 'Sequence.pm';



$db[0] = "/home/olly/COMPARATIVE_YEAST/ORFS/BLAST_REPORTS/BINDINGMAP/DATABASES/UTR5/utr5_sc_1000.fasta";
$db[3] =  "/home/olly/YEAST_GENOMES/utr5_1000_mikatae_mit.fasta";
$db[2] =  "/home/olly/YEAST_GENOMES/utr5_1000_bayanus_mit.fasta";
$db[1] =  "/home/olly/YEAST_GENOMES/utr5_1000_paradoxus_mit.fasta";

$d = ($ARGV[2]?$ARGV[2]:0);
#print "d=$d\n";
my $o_seq = Sequence->new;
$o_seq->setBlastPath("/home/olly/COMPARATIVE_YEAST/ORFS/BLAST_REPORTS/BINDINGMAP/PROGRAMS/BLAST");
$o_seq->setBlastDB($db[$d]);

my $s1 = $o_seq->getSequenceFromBlastDB($ARGV[1], 0, 0);



if ($s1) {
    open OUT, ">tmpo";
    print OUT ">$ARGV[1]\n$s1\n\n";
    close OUT;
    
} else {
    exit;
}


$s = ScanACE->new;

#$s->readMotif("../MOTIFS/RAP1.motif");
$s->setNbMotifs(100);
#$s->setStdDev(3.0);
$s->setGC(0.33);
#$s->setRun(0);
#$s->setThresholdScore(10.0);

$s->runScanACE("tmpo", $ARGV[0]);


#my $cnt = $s->getNbMotifsInFile($ARGV[0]);

my $a_ref = $s->getSites($i);

$tot = 0.0;    
foreach my $a (@$a_ref) {
    
    print "$ARGV[0]\t$a->[0]\t$a->[3]\t$a->[4]\t$a->[1]\t$a->[2]\n";
    
    $tot += $a->[4];
}


    #print "tot=$tot\n";
