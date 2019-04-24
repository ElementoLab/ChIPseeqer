use lib qw(/home/elemento/PERL_MODULES);

use Fasta;
use Sets;
use MyBlast;
use strict;



my $fa = Fasta->new;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");

$mb->setDatabaseDatabase($ARGV[1]);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("1e-3");
$mb->setMegablast(1);

    
#
# go thru all the sequences, align them one by one
#

my $file = Sets::getTempFile("toto");


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    $fa->writeSeq($file, $n, $s);
    $mb->setQueryFile($file);
    
    my $a_ref = $mb->blastallMultiple;

    foreach my $hit (@$a_ref) { 
	
	my $hsps = $hit->{"HSPS"};
	my $h    = $hit->{"HIT_NAME"};
	foreach my $r (@$hsps) {
	    
	    my $l  = $r->{"ALIGNLEN"};
	    # my $i  = $r->{"IDENTITY"};
	    
	    my $s  = $r->{"QFROM"};
	    my $e  = $r->{"QTO"};

	    my $ev = $r->{"EVALUE"};
	    my $st = $r->{"DFRAME"};

	    my $i  = Sets::getSequencesIdentity($r->{DSEQ}, $r->{QSEQ});

	    if (($l >= 100) && ($i >= 0.8)) {
	      print "$n\t$l\t$i\t$e\t$s\n";
	    }
	}
    }

}

unlink $file;

