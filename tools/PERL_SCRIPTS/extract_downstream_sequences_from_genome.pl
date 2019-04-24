BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Table;
use Sequence;
use Sets;
use strict;

if (scalar(@ARGV) == 0) {
	print "Usage : extract_downstream_sequences_from_genome.pl annotation genome lengthD lengthU minlen [ 3'UTR lengths ]\n";
	exit(0);
}
my $t = Table->new;
$t->loadFile($ARGV[0]);


my $s = Sequence->new;
$s->setVerbose(1);
$s->setBlastDB($ARGV[1]);
my $lenD   = $ARGV[2];
my $lenU   = $ARGV[3];
my $minlen = $ARGV[4];



my $a_ref = $t->getArrayOfHashes( ("ORF", "SCAFFOLD", "START_P", "END_P", "STRAND", "START_T", "END_T") );


#
#  load 3'UTR lengths to add if available
#
my $h_ref_len = undef;

if (defined($ARGV[5])) {
    my %HH = ();
    $h_ref_len = \%HH;
    
    $t->loadFile($ARGV[5]);
    $h_ref_len = $t->getIndexKV(0,1);
}

 

foreach my $r (@$a_ref) {

    if ($r->{"END_T"} < $r->{"START_T"}) {
	my $tt = $r->{"START_T"};
	$r->{"START_T"} = $r->{"END_T"};
	$r->{"END_T"} = $tt;
    }

    if ($r->{"END_P"} < $r->{"START_P"}) {
	my $tt = $r->{"START_P"};
	$r->{"START_P"} = $r->{"END_P"};
	$r->{"END_P"} = $tt;
    }

    
    #
    #  case where only the protein sequences are defined
    #
    if (!defined($r->{"END_T"})) {

	$r->{"START_T"} = $r->{"START_P"};
	$r->{"END_T"}   = $r->{"END_P"};

	#
	#  case where we want to add a 3'UTR length
	#

	if (defined($h_ref_len) && (defined($h_ref_len->{ $r->{ORF} }))) {
	    
	    if ($r->{STRAND} < 0) {
		$r->{START_T} = Sets::max(0, $r->{START_T} - $h_ref_len->{ $r->{ORF} });
	    } else {
		$r->{END_T}   = $r->{END_T}  + $h_ref_len->{ $r->{ORF} };
	    }	    
	}
	
    }
    

    my $start     = undef;
    my $end       = undef;

    
    if ($r->{STRAND} < 0) {
	
	$start = Sets::max(0, $r->{START_T} - $lenD);
	$end   = $r->{START_T} + $lenU;
    
    } else {
	
	$start = Sets::max(0, $r->{END_T} - $lenU);
	$end   = $r->{END_T} + $lenD;
	
    }

    next if ( ($end - $start) <= $minlen);
     
    my $seq = $s->getSequenceFromBlastDB($r->{SCAFFOLD}, $start, $end);
    
    if ($r->{STRAND} < 0) {
	$seq = Sets::getComplement($seq);
    }


    if ($seq && (length($seq) > $minlen)) {
        print ">$r->{ORF}|$r->{SCAFFOLD}|$start|$end|$r->{STRAND}\n$seq\n\n";
    }
    
    



}
