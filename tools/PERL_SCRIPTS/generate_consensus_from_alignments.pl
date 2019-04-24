use lib qw(/home/olly/PERL_MODULES);
use Sets;
use  ClustalW;

my $a_ref_files = Sets::readSet($ARGV[0]);

foreach my $f (@$a_ref_files) {
    my $cl = ClustalW->new;
    $cl->setFile($f);
    
    my $a_ref_aln = $cl->getSeqArray;
    
    if (scalar(@$a_ref_aln) >= $ARGV[1]) {
	system("perl ~/PERL_MODULES/SCRIPTS/aln2consensus.pl $f > $ARGV[2]/$f");
    }
}

