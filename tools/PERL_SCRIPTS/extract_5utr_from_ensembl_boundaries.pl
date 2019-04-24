use lib qw(/home/elemento/PERL_MODULES);

use Table;
use Sequence;
use Sets;
use strict;

if (scalar(@ARGV) == 0) {
	print "Usage : extract_upstream_sequences_from_blast_predictions.pl blast_matches.txt DB2 LEN [du c? ] [lens maxlen]\n";
}
my $t = Table->new;
$t->loadFile($ARGV[0]);


my $s = Sequence->new;
#$s->setBlastPath("/home/olly/COMPARATIVE_YEAST/ORFS/BLAST_REPORTS/BINDINGMAP/PROGRAMS/BLAST");
$s->setBlastDB($ARGV[1]);
#$s->setVerbose(1);
my $len = $ARGV[2];





my $a_ref = $t->getArrayOfHashes( ("ORF", "SCAFFOLD", "START_P", "END_P", "STRAND", "START_T", "END_T") );

my %HH = ();
my $h_ref_len = \%HH;

if (defined($ARGV[3])) {
 $t->loadFile($ARGV[3]);
 $h_ref_len = $t->getIndexKV(0,1);

} 

foreach my $r (@$a_ref) {

    my $start = undef;
    my $end   = undef;
 
    if ($r->{"START_P"} > $r->{"END_P"}) {
      my $tmp = $r->{"START_P"};
      $r->{"START_P"} = $r->{"END_P"};
      $r->{"END_P"}   = $tmp;
    }

    if (!defined($r->{"START_T"})) {
      $r->{"START_T"} = $r->{"START_P"};
    }

    if (!defined($r->{"END_T"})) {
      $r->{"END_T"} = $r->{"END_P"};
    }
	
    #print join("\t", values(%$r)); print "\n";
    #print join("\t", keys  (%$r)); print "\n";
    
    # what do we want ?
    if ($r->{STRAND} < 0) {
	
	$start = $r->{"END_P"};
	$end   = $r->{"END_T"};

	if (defined($h_ref_len->{ $r->{ORF} }) && ($start == $end)) {
	  $end += $h_ref_len->{ $r->{ORF} };
	} elsif (defined($len) && ($start == $end)) {
	  $end += $len;
	}

	if (abs($end - $start) > 5000) {
	  $end   = $r->{"END_P"} + 5000;
	}

    } else {
	
	$start = $r->{"START_T"};
	$end   = $r->{"START_P"}; 
	
	if (defined($h_ref_len->{ $r->{ORF} }) && ($start == $end)) {
	  $start -= $h_ref_len->{ $r->{ORF} };
	} elsif (defined($len) && ($start == $end)) {
	  $start -= $len; 
	}

	if (abs($end - $start) > 5000) {
	  $start   = $r->{"START_P"} - 5000;
	}

	if ($start < 0) {
	    $start = 1;
	}

	
	
    }
   
    #print "$r->{ORF}, $start, $end\n";
     
    next if (abs($start - $end) <= 5);

    my $seq = $s->getSequenceFromBlastDB($r->{SCAFFOLD}, $start, $end);
    
    if ($r->{STRAND} < 0) {
	$seq = Sets::getComplement($seq);
    }

    if ($seq && (length($seq) > 1)) {
        print ">$r->{ORF}\n$seq\n\n";
    }
    
    

}
