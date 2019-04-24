use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;
use Sequence;
use Fasta;
use strict;
use MyBlast;


my $mb = MyBlast->new;
$mb->setBlastProgram("tblastn");
$mb->setDatabaseDatabase($ARGV[1]);
my $tmpfile1 = Sets::getTempFile("/tmp/blast.1");
my $tmpfile2 = Sets::getTempFile("/tmp/blast.2");


$mb->setEvalueThreshold("0.1");
$mb->setNbProcessors(2);
$mb->setVerbose(0);


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
#  load the dmel genes
#

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
    
    my $k1 = $r->[0];

    #
    # traverse all the orthologs
    #

    my $cnt = 0;
    foreach my $o (@a_orthologs) {
	$cnt++ if (defined($o->{ORTHOLOGS}->{ $k1 }));
    }


    next if ($cnt != scalar(@a_orthologs));
    
    
    #
    #  GET the aa sequence
    #
    $se->setBlastDB($ARGV[1]);
    my $seqaa = $se->getSequenceFromBlastDB($k1, 0, 0);

    $fa->writeSeq($tmpfile1, $k1, $seqaa);

    #
    # traverse all the orthologs again ..
    #
    
    my $flag_zero_length = 0;
    my $txt = "";
    foreach my $o (@a_orthologs) {
	
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

	my $seq_orth = $se->getSequenceFromBlastDB($ch, $st, $en);
	
	
	#
	#  blast the dmel aa sequences againt the nt sequences
	#
	$fa->writeSeq($tmpfile2, $k1, $seq_orth);

	$mb->setBlastProgram("tblastn");

	$mb->setQueryDatabase($tmpfile1);
	$mb->setDatabaseFile($tmpfile2);
	
	my $a_ref_hsp = $mb->blastallUnique;
	

	#
	#  start collecting the HSP
	#
	my @a_pieces = ();
	
	foreach my $hsp (@$a_ref_hsp) {
    
	    # the current HSP is compatible if it does not overlap too much with any current piece
	    my $overlap = 0;
	    foreach my $p (@a_pieces) {
		$overlap = Sets::getSequencesOverlap($hsp->{QFROM}, $hsp->{QTO}, $p->{QFROM}, $p->{QTO});
		last if ($overlap > 30);
	    }
	    
	    if ($overlap < 30) {
		push @a_pieces, $hsp;
	    }
	}
	@a_pieces = sort { $a->{QFROM} <=> $b->{QFROM} } @a_pieces;
    
	
	#
	#  get the contig / start / end 
	#
	my $d_start = ($a_pieces[0         ]->{DFRAME}>0?$a_pieces[0         ]->{DFROM}:$a_pieces[0         ]->{DTO});
	my $d_end   = ($a_pieces[$#a_pieces]->{DFRAME}>0?$a_pieces[$#a_pieces]->{DTO}:$a_pieces[$#a_pieces]->{DFROM});

	#
	#  get the protein sequence
	#
	
	my $frameN = 0;
	my $frameP = 0;
	
	my $orth_protein = "";
	
	my $i = 0;
	foreach my $p (@a_pieces) {
	    $orth_protein .= $p->{DSEQ};
	    if ($p->{DFRAME} < 0) {
		$frameN += 1;
	    } else {
		$frameP += 1;
	    }
	    $i++;
	}
	
	next if (($frameP >  0) && ($frameN >  0)); 
	next if (($frameP == 0) && ($frameN == 0));
	
	my $frame = undef;
	if ($frameP > $frameN) {
	    $frame =  1;
	} else {
	    $frame = -1;
	}
	
	
	
	$orth_protein =~ s/\-//g;
	my $length = length($orth_protein);
	
	if ($length == 0) {
	    $flag_zero_length = 1;
	    last;
	}

	if ($orth_protein =~ /\*/) {
	    $flag_zero_length = 1;
	    last;
	}

	$txt .= ">$o->{NAME}\n$orth_protein\n\n";
	
    } 

    next if ($flag_zero_length == 1);

    $txt = ">dmel\n$seqaa\n\n" . $txt;

    open OUT, ">PROTEINS/$k1.seq";
    print OUT $txt;
    close OUT;
    
}
