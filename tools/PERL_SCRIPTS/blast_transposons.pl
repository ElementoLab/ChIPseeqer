use lib qw(/home/olly/PERL_MODULES);

use Fasta;
use Sets;
use MyBlast;
use strict;



my $fa = Fasta->new;
$fa->setFile($ARGV[0]);


my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");

$mb->setDatabaseDatabase($ARGV[1]);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("1e-50");

#$mb->setMismatchWeight(-6);
#$mb->setMatchWeight(5);
#$mb->setGapOpening(20);
#$mb->setGapExtension(4);
#$mb->setWordLength(7);
#$mb->setQueryStrand(1);
    
    
#
# go thru all the sequences, align them one by one
#

my $file = Sets::getTempFile("toto");
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    $fa->writeSeq($file, $n, $s);
    $mb->setQueryFile($file);
    


    
    my $a_ref = $mb->blastallMultiple;

    #print "got " . scalar(@$a_ref) . " hits\n";

    foreach my $r (@{ $a_ref->[0]->{"Hit_hsps"} }) {
	
	my $s1 = $r->{"Hsp_hit-from"};
	my $e1 = $r->{"Hsp_hit-to"};
	my $ev = $r->{"Hsp_evalue"};
	my $st = $r->{"Hsp_hit-frame"};

	
	$n =~ s/\#/\t/;

	print "$ARGV[1]\t$n\t$s1\t$e1\t$ev\t$st\n";
    }

    #<STDIN>;
}

unlink $file;
