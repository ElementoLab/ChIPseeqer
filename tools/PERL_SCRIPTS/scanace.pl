use lib qw(/home/olly/PERL_MODULES);
require "ScanACE.pm";

$s = ScanACE->new;

#$s->readMotif("../MOTIFS/RAP1.motif");
#$s->setNbMotifs(5000);
$s->setStdDev(3.0);
$s->setGC(0.33);
#$s->setRun(0);
#$s->setThresholdScore(10.0);

$s->runScanACE($ARGV[1], $ARGV[0]);


my $cnt = $s->getNbMotifsInFile($ARGV[0]);

    my $a_ref = $s->getSites($i);
    
    
    foreach my $a (@$a_ref) {
	
	print "$ARGV[0]\t$a->[0]\t$a->[3]\t$a->[4]\t$a->[1]\t$a->[2]\n";
	
    }

