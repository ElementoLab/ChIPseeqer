use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;
use Sequence;
use Fasta;
use strict;


my $fa = Fasta->new;

my $se = Sequence->new;


my $ta = Table->new;

#
#  load the orthologs
#
$ta->loadFile("ortholog_table.txt");
my $a_ref_orth_info = $ta->getArray();
my @a_orthologs = ();
foreach my $r (@$a_ref_orth_info) {
    
    $ta->loadFile($r->[2]);
    
    my $h_ref_o = $ta->getIndex(0);

    my %h_tmp = (ORTHOLOGS => $h_ref_o, FILE => $r->[1], NAME => $r->[0]);
    
    push @a_orthologs, \%h_tmp;

}




#
#  load the exons
#

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %CHR = ();
my %EXONS = ();
my %BOUNDARIES = ();
my %STRAND = ();
foreach my $r (@$a_ref) {
    my @a_exon = ($r->[2], $r->[3]);
    push @{ $EXONS{ $r->[1] } }, \@a_exon;
    $CHR{ $r->[1] } = $r->[0]; 
    $BOUNDARIES{ $r->[1] }->[0] = (defined($BOUNDARIES{ $r->[1] }->[0])?Sets::min($BOUNDARIES{ $r->[1] }->[0], $r->[2]):$r->[2]);
    $BOUNDARIES{ $r->[1] }->[1] = (defined($BOUNDARIES{ $r->[1] }->[1])?Sets::max($BOUNDARIES{ $r->[1] }->[1], $r->[3]):$r->[3]);
    $STRAND{ $r->[1] } = $r->[4];
}

my @keys = keys(%EXONS);


my $w    = 5000;
foreach my $k1 (@keys) {
    
    
    
    next if ($k1 ne $ARGV[1]);

    #print "--> $k1 $STRAND{$k1}\n";

    # 
    #  get the coding and surrounding regions
    #
    my $start = $BOUNDARIES{ $k1 }->[0] - $w;
    my $end   = $BOUNDARIES{ $k1 }->[1] + $w;
    $se->setBlastDB("dmel.fasta");
    my $seq = $se->getSequenceFromBlastDB($CHR{ $k1}, $start, $end);

   
    #
    # mask the exons
    #
    my @a_exons = ();
    foreach my $e (@{ $EXONS{ $k1 } }) {
	my $start_exon = $e->[0] -  $BOUNDARIES{ $k1 }->[0] + $w;
	my $end_exon   = $e->[1] -  $BOUNDARIES{ $k1 }->[0] + $w;
	my @a          = ($start_exon, $end_exon);
	push @a_exons, \@a;
    }    
    my $seq_masked = Sets::maskExons($seq, \@a_exons, 'N');

    
    # get all the overlaping guys, within a +- 5000 nt window
    foreach my $k2 (@keys) {
	
	next if ($k1 eq $k2);
	next if ($CHR{$k1} ne $CHR{$k2});
	next if (!Sets::sequencesOverlap($BOUNDARIES{ $k1 }->[0]-$w, $BOUNDARIES{ $k1 }->[1]+$w, $BOUNDARIES{ $k2 }->[0], $BOUNDARIES{ $k2 }->[1]) );
	
	# get only the overlapping exons
	my @a_exons = ();
	foreach my $e (@{ $EXONS{ $k2 } }) {
	    next if (!Sets::sequencesOverlap($BOUNDARIES{ $k1 }->[0]-$w, $BOUNDARIES{ $k1 }->[1]+$w, $e->[0], $e->[1]) );
	    my $start_exon = Sets::max(1, $e->[0] -  $BOUNDARIES{ $k1 }->[0] + $w);
	    my $end_exon   = Sets::min(length($seq_masked), $e->[1] -  $BOUNDARIES{ $k1 }->[0] + $w);

	    

	    #my $end_exon   = $e->[1] -  $BOUNDARIES{ $k1 }->[0] + $w;
	    my @a          = ($start_exon, $end_exon);
	    push @a_exons, \@a;
	}    

	$seq_masked = Sets::maskExons($seq_masked, \@a_exons, 'X');
		
    }

    my $fa = Fasta->new;

    my $tmpfile = "OUT/$k1" . "_dmel.seq";
    $fa->writeSeq($tmpfile, "$k1" . "_dmel", $seq_masked);

    #
    # get the same regions in orthologous genomes 
    #

    open OUT, ">OUT/$k1" . "_dorth.seq";
    foreach my $o (@a_orthologs) {

	next if (!defined($o->{ORTHOLOGS}->{ $k1 }));
	
	my $st = $o->{ORTHOLOGS}->{ $k1 }->[ 2 ];
	my $en = $o->{ORTHOLOGS}->{ $k1 }->[ 3 ];

	if ($en < $st) {
	    my $tm = $en;
	    $en = $st;
	    $st = $tm;
	}

	my $ch = $o->{ORTHOLOGS}->{ $k1 }->[ 1 ];

	
	if ($o->{NAME} eq "dana") {
	    $ch = "C" . $ch;
	}

	$se->setBlastDB($o->{FILE});
	#$se->setVerbose(1);

	my $seq_orth = undef;
	if ($o->{ORTHOLOGS}->{ $k1 }->[ 4 ] == $STRAND{ $k1 }) { 
	    $seq_orth = $se->getSequenceFromBlastDB($ch, Sets::max(0,$st-$w), $en+$w);
	} else {

	    $seq_orth = $se->getSequenceFromBlastDB($ch, Sets::max(0,$st-$w), $en+$w);
	    $seq_orth = Sets::getComplement($seq_orth);
	}
	
	print OUT ">$k1 $o->{FILE}\n$seq_orth\n";

    } 
    close OUT;



    # trf ..


    # local alignment
    system("/usr/bin/perl /home/olly/PERL_MODULES/SCRIPTS/do_lalign.pl OUT/$k1" . "_dmel.seq OUT/$k1" . "_dorth.seq OUT_RES/$k1.txt");

    
    
    
    
}
