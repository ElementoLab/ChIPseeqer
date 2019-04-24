use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;
use Sequence;
use Fasta;
use strict;
use Repeats;
#use Data::Dumper;

my $fa = Fasta->new;

my $se = Sequence->new;


my $ta = Table->new;


#
#  load the exons
#

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %CHR            = ();
my %EXONS          = ();
my %BOUNDARIES     = ();
my %BOUNDARIES_CDS = ();
my %STRAND         = ();

foreach my $r (@$a_ref) {
    my @a_exon = ($r->[2], $r->[3]);
    push @{ $EXONS{ $r->[1] } }, \@a_exon;
    $CHR{ $r->[1] }             = $r->[0]; 
    $BOUNDARIES{ $r->[1] }->[0] = (defined($BOUNDARIES{ $r->[1] }->[0])?Sets::min($BOUNDARIES{ $r->[1] }->[0], $r->[2]):$r->[2]);
    $BOUNDARIES{ $r->[1] }->[1] = (defined($BOUNDARIES{ $r->[1] }->[1])?Sets::max($BOUNDARIES{ $r->[1] }->[1], $r->[3]):$r->[3]);
    $STRAND{ $r->[1] }          = $r->[4];

}


my @keys = keys(%EXONS);


my $wmax    = 20000;

my $cnt     = 0;
foreach my $k1 (@keys) {
    
    $cnt ++;

    #next if ($cnt < 1056);
    
    #next if ($k1 ne $ARGV[1]);

    #print "--> $k1 $STRAND{$k1}\n";

    # 
    #  get the coding and surrounding regions
    #

    my $w = ($BOUNDARIES{ $k1 }->[0] - $wmax < 0?$BOUNDARIES{ $k1 }->[0]-1:$wmax);


    my $start = Sets::max(0, $BOUNDARIES{ $k1 }->[0] - $w);
    my $end   = $BOUNDARIES{ $k1 }->[1] + $w;  # in principle should use min here
    $se->setBlastDB($ARGV[1]);

    #print "getting $CHR{ $k1}, $start, $end\n";
    
    my $seq_masked = $se->getSequenceFromBlastDB($CHR{ $k1}, $start, $end);
    
    next if (!$seq_masked);

    my $l = length($seq_masked);
    
    #
    # mask the OTHER exons
    #
    foreach my $k2 (@keys) {
	
	next if ($k1 eq $k2);
	next if ($CHR{$k1} ne $CHR{$k2});
	next if (!Sets::sequencesOverlap($BOUNDARIES{ $k1 }->[0]-$w, $BOUNDARIES{ $k1 }->[1]+$w, $BOUNDARIES{ $k2 }->[0], $BOUNDARIES{ $k2 }->[1]) );
	
	# get only the overlapping exons
	my @a_exons = ();
	foreach my $e (@{ $EXONS{ $k2 } }) {
	    next if (!Sets::sequencesOverlap($BOUNDARIES{ $k1 }->[0]-$w, $BOUNDARIES{ $k1 }->[1]+$w, $e->[0], $e->[1]) );

	    # get the coordinates of the exons, in the $w + coding + $w reference frame
	    my $start_exon = Sets::max(1, $e->[0] -  $BOUNDARIES{ $k1 }->[0] + $w);
	    my $end_exon   = Sets::min(length($seq_masked), $e->[1] -  $BOUNDARIES{ $k1 }->[0] + $w);

	    
	    #print "masking $start_exon to $end_exon\n";

	    my @a          = ($start_exon, $end_exon);

	    #print EXO "$k2\t$start_exon\t$end_exon\n";
	    
	    push @a_exons, \@a;
	}    

	$seq_masked = Sets::maskExons($seq_masked, \@a_exons, 'N');
	
	#print "masking " . scalar(@a_exons) . " exons from $k2\n";
 
    }

    #
    #  mask the real EXONS
    # 
    my @a_exons = ();
    foreach my $e (@{ $EXONS{ $k1 } }) {
	my $start_exon = $e->[0] -  $BOUNDARIES{ $k1 }->[0] + $w;
	my $end_exon   = $e->[1] -  $BOUNDARIES{ $k1 }->[0] + $w;
	
	#print "masking $start_exon to $end_exon in seq of length $l\n";

	
	my @a          = ($start_exon, $end_exon);
	push @a_exons, \@a;
    }    
    my $seq_masked = Sets::maskExons($seq_masked, \@a_exons, 'X');

    
    #
    # get all the overlaping guys, within a +- $w nt window 
    #
  
    #print "$l_utr5 $l_utr3\n";

    #system("mkdir RESULTS/$k1") if (! -e "RESULTS/$k1");

    
    #open EXO, ">RESULTS/$k1/dmel.other_exons";

    
    #close EXO;
    
    #print $seq_masked;

    #
    #  mask the tandem repeats with period <= 3
    #
    #my $tr = Repeats->new;
    #$tr->setSeq($seq_masked);
    #$tr->setSpecies("primates");
    #$tr->setProgram('REPEATMASKER');
    #$tr->process();
    #my $a_ref_tandems = $tr->getResults;
    #$tr->dispose;

    #$seq_masked = Sets::maskExons($seq_masked, $a_ref_tandems, 'X');

    #print "masking " . scalar(@$a_ref_tandems) . " repeats\n";

    if ($STRAND{$k1} == -1) {
	$seq_masked = getSpeComplement($seq_masked);
    }


    print ">$k1\n$seq_masked\n\n";

    
}


sub getSpeComplement {
    
    my ($str) = @_;
    
    my $l = length($str);
    
    my @s = split //, $str;

    my $c = "";
    for (my $i=$l-1; $i>=0; $i--) {
	
	my $d = "";
	
	if ($s[$i] eq 'A') {
	    $d = 'T';
	} elsif ($s[$i] eq 'T') {
	    $d = 'A';
	} elsif ($s[$i] eq 'G') {
	    $d = 'C';
	} elsif ($s[$i] eq 'C') {
	    $d = 'G';
	} elsif ($s[$i] eq '-') {
	    $d = '-';
	} else {
	    $d = $s[$i];
	}

	$c .= $d;
    }
    
    return $c;
   
}
