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
#  load the orthologs
#
$ta->loadFile("ortholog_table.txt");
my $a_ref_orth_info = $ta->getArray();
my @a_orthologs = ();
my $nb = 1;
foreach my $r (@$a_ref_orth_info) {
    
    $ta->loadFile("DATA/$r->[2]");
    
    my $h_ref_o = $ta->getIndex(0);

    my %h_tmp = (ORTHOLOGS => $h_ref_o, FILE => $r->[1], NAME => $r->[0], NUMBER => $nb);
    
    push @a_orthologs, \%h_tmp;

    $nb ++;
}




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

    # no need to do that if this exon is not part of the CDS
    next if (($r->[5] eq "") && ($r->[6] eq ""));

    #print Dumper($r);
    
    $BOUNDARIES_CDS{ $r->[1] }->[0] = (defined($BOUNDARIES_CDS{ $r->[1] }->[0])?Sets::min($BOUNDARIES_CDS{ $r->[1] }->[0], $r->[5]):$r->[5]);
    $BOUNDARIES_CDS{ $r->[1] }->[1] = (defined($BOUNDARIES_CDS{ $r->[1] }->[1])?Sets::max($BOUNDARIES_CDS{ $r->[1] }->[1], $r->[6]):$r->[6]);
}


#print Dumper(%BOUNDARIES_CDS);

my @keys = keys(%EXONS);


my $w    = 10000;
foreach my $k1 (@keys) {
    
    
    
    next if ($k1 ne $ARGV[1]);

    print "--> $k1 $STRAND{$k1}\n";

    # 
    #  get the coding and surrounding regions
    #
    my $start = Sets::max(0, $BOUNDARIES{ $k1 }->[0] - $w);
    my $end   = $BOUNDARIES{ $k1 }->[1] + $w;  # in principle should use min here
    $se->setBlastDB("DATA/dmel.fasta");
    #$se->setVerbose(1);
    my $seq = $se->getSequenceFromBlastDB($CHR{ $k1}, $start, $end);

   
    #
    # mask the exons
    #

    my $l = length($seq);
    
    print "sequence to align is $l bp\n";
    
    my @a_exons = ();
    foreach my $e (@{ $EXONS{ $k1 } }) {
	my $start_exon = $e->[0] -  $BOUNDARIES{ $k1 }->[0] + $w;
	my $end_exon   = $e->[1] -  $BOUNDARIES{ $k1 }->[0] + $w;

	#print "masking $start_exon to $end_exon\n";

	
	my @a          = ($start_exon, $end_exon);
	push @a_exons, \@a;
    }    
    my $seq_masked = Sets::maskExons($seq, \@a_exons, 'N');

    
    #
    # get all the overlaping guys, within a +- $w nt window 
    #
  
    #print "$l_utr5 $l_utr3\n";

    system("mkdir RESULTS/$k1") if (! -e "RESULTS/$k1");

    
    open EXO, ">RESULTS/$k1/dmel.other_exons";

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

	    print EXO "$k2\t$start_exon\t$end_exon\n";
	    
	    push @a_exons, \@a;
	}    

	$seq_masked = Sets::maskExons($seq_masked, \@a_exons, 'X');
	
	print "masking " . scalar(@a_exons) . " exons from $k2\n";
 
    }
    close EXO;
    

    #
    #  mask the tandem repeats with period <= 3
    #
    my $tr = Repeats->new;
    $tr->setSeq($seq_masked);
    $tr->setSpecies("drosophila");
    $tr->setProgram('REPEATMASKER');
    $tr->process();
    my $a_ref_tandems = $tr->getResults;
    $tr->dispose;

    open TAN, ">RESULTS/$k1/dmel.repeats";
    foreach my $r (@$a_ref_tandems) {
	print TAN "$r->[0]\t$r->[1]\t$r->[2]\n";
    }
    close TAN;
    
    #print "BEFORE " . length($seq_masked);

    $seq_masked = Sets::maskExons($seq_masked, $a_ref_tandems, 'X');

    print "masking " . scalar(@$a_ref_tandems) . " repeats\n";

    #print "AFTER " . length($seq_masked);
    #<STDIN>;
    #
    #  mask the repeats using repeatmasker
    #
    

    my $fa = Fasta->new;

    my $tmpfile = "RESULTS/$k1/dmel.seq";
    $fa->writeSeq($tmpfile, "$k1" . "_dmel|$CHR{$k1}|$start|$STRAND{$k1}", $seq_masked);


    my $tmpfile = "RESULTS/$k1/dmel.seq.init";
    $fa->writeSeq($tmpfile, "$k1" . "_dmel|$CHR{$k1}|$start|$STRAND{$k1}", $seq);

    #
    # get the same regions in orthologous genomes 
    #
    
    my $l_utr5 = abs($BOUNDARIES{ $k1 }->[0] - $BOUNDARIES_CDS{ $k1 }->[0]);
    my $l_utr3 = abs($BOUNDARIES{ $k1 }->[1] - $BOUNDARIES_CDS{ $k1 }->[1]);

    print "5'utr = $l_utr5 bp, 3'utr = $l_utr3\n";

    open OUT, ">RESULTS/$k1/dort.seq";
    foreach my $o (@a_orthologs) {

	next if (!defined($o->{ORTHOLOGS}->{ $k1 }));

	#print "Ortholog exists in $o->{NAME}\n";
	
	my $st = $o->{ORTHOLOGS}->{ $k1 }->[ 2 ];
	my $en = $o->{ORTHOLOGS}->{ $k1 }->[ 3 ];

	if ($en < $st) {
	    my $tm = $en;
	    $en = $st;
	    $st = $tm;
	}

	

	my $ch = $o->{ORTHOLOGS}->{ $k1 }->[ 1 ];
	my $len = $en-$st;
	print "Ortholog exists in $o->{NAME} ($len bp, chr is $ch)\n";
	
	#if ($o->{NAME} eq "dana") {
	#    $ch = "C" . $ch;
	#}

	$se->setBlastDB("DATA/$o->{FILE}");
	#$se->setVerbose(1);

	my $seq_orth = undef;
	if ($o->{ORTHOLOGS}->{ $k1 }->[ 4 ] == $STRAND{ $k1 }) { 
	    $seq_orth = $se->getSequenceFromBlastDB($ch, Sets::max(0,$st-$w-$l_utr5), $en+$w+$l_utr3);
	} else {

	    #  inverse strand, need to invert also the utrs
	    $seq_orth = $se->getSequenceFromBlastDB($ch, Sets::max(0,$st-$w-$l_utr3), $en+$w+$l_utr5);
	    $seq_orth = Sets::getComplement($seq_orth);
	}
	
	print OUT ">$k1|$o->{NAME}|$o->{NUMBER}\n$seq_orth\n";

	#print "out ..\n";
    } 
    close OUT;

    print "align now ..\n";
    # local alignment
    system("/usr/bin/perl /home/olly/PERL_MODULES/SCRIPTS/do_lalign_tmp.pl RESULTS/$k1/dmel.seq RESULTS/$k1/dort.seq");

    system("/home/olly/PROGRAMS/ENHANCERS/calc_conservation_statistics -consfile RESULTS/$k1/dmel.seq.cons -w 200 -minlen 200  -nbstd 3.0 -outfile RESULTS/$k1/dmel.seq.scores > RESULTS/$k1/dmel.seq.enh");
    
    
    
}
