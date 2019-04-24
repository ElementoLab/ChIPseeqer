use lib qw(/home/olly/PERL_MODULES);

use MotifLibrary;

my $mo = MotifLibrary->new;
$mo->setMaxNbMotifs(10);
$mo->loadAlignACEOutputFile($ARGV[0]);

$mo->setMaxNbMotifs(20);
$mo->loadAlignACEOutputFile($ARGV[1]);
 

$mo->sortByScore();
#$mo->printInfo();

my $a_ref = $mo->getMotifs();

my $nbm = 1;
foreach my $m (@$a_ref) {
    
    #$m->printAlignACE;
    my $size    = $m->getSize();
    my $cnt     = $m->getNbColumnsWithMaxNbSymbols(2);
    my $cnt_seq = $m->getNumberDistinctSequences(0);
    #print "$cnt/$size columns with < 2 or less nt, $cnt_seq for droso\n";

    if (($cnt/$size > 0.6) && ($cnt_seq >= 4)) {
	print "Motif $nbm\n";
	$m->printAlignACE;
	print "\n";
	#print "$cnt/$size columns with < 2 or less nt, $cnt_seq for droso\n";
	$nbm++;
    }
}
