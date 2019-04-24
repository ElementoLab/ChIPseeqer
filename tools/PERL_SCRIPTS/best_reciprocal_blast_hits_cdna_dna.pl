#
# input : set of proteins vs set of proteins
#
use lib qw(/home/olly/PERL_MODULES);
use MyBlast;
use Fasta;
use Sets;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");
$mb->setDatabaseDatabase($ARGV[1]);
my $tmpfile = "/tmp/blast.1";
$mb->setQueryDatabase($tmpfile);
$mb->setEvalueThreshold("1e-10");
#$mb->setMismatchWeight(-);
#$mb->setVerbose(1);

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $qlen = undef;

while (my $a_ref = $fa->nextSeq) {
    
    my ($name, $seq) = @$a_ref;

    my $qlen = length($seq);

    #print "seq1=$seq\n";

    my @a = split / /, $name;

    # create a query file
    $fa->writeSeq($tmpfile, $name, $seq);

    # do the blast
    #$mb->setVerbose(1);
    $mb->setBlastProgram("blastn");
    $mb->setDatabaseDatabase($ARGV[1]);
    $mb->setQueryDatabase($tmpfile);

    my $a_ref = $mb->blastallUnique;
    

    next if (scalar(@$a_ref) == 0);

    # put the pieces back together, from high scoring to low scoring
    
    my @a_pieces = ();
    
    foreach my $r (@$a_ref) {
    
	#print "$r->{EVALUE}\t$r->{QFROM}\t$r->{QTO}\t$r->{DFRAME}\n$r->{QSEQ}\n";  

	# the current HSP is compatible if it does not overlap too much with any current piece
	
	my $overlap = 0;
	foreach my $p (@a_pieces) {
	    
	    $overlap = Sets::getSequencesOverlap($r->{QFROM}, $r->{QTO}, $p->{QFROM}, $p->{QTO});

	    #print "overlap between $r->{QFROM}, $r->{QTO} and  $p->{QFROM}, $p->{QTO} is $overlap\n";

	    last if ($overlap > 3);
	}

	if ($overlap <= 3) {
	    push @a_pieces, $r;
	}
	
	#substr($seq, $r->{QFROM} - 1, $r->{QTO} - $r->{QFROM} + 1)  = 'N' x ( $r->{QTO} - $r->{QFROM} + 1);
	
	#print "\n";
    }

    
    @a_pieces = sort { $a->{QFROM} <=> $b->{QFROM} } @a_pieces;
    

    #
    #  get the contig / start / end 
    #
    my $d_start = ($a_pieces[0         ]->{DFRAME}>0?$a_pieces[0         ]->{DFROM}:$a_pieces[0         ]->{DTO});
    my $d_end   = ($a_pieces[$#a_pieces]->{DFRAME}>0?$a_pieces[$#a_pieces]->{DTO}:$a_pieces[$#a_pieces]->{DFROM});
    my $d_id    = $mb->getUniqueHitName();


    #print $a_pieces[$#a_pieces]->{DFRAME}; <STDIN>;

    #
    #  get the protein sequence
    #

    my $frameN = 0;
    my $frameP = 0;

    my $orth_seq = "";

    my $i = 0;
    foreach my $r (@a_pieces) {
	
	#print "$r->{EVALUE}\t$r->{QFROM}\t$r->{QTO}\t$r->{DFROM}\t$r->{DTO}\t$r->{DFRAME}\n$r->{QSEQ}\n$r->{DSEQ}\n";  
	#print "\n";
	
	#$length += ($r->{QTO} - $r->{QFROM} + 1);
 
	$orth_seq .= $r->{DSEQ};
	
	if ($r->{DFRAME} < 0) {
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
    
    

    $orth_seq =~ s/\-//g;
    my $length = length($orth_seq);

    
    #print "Total length aligned = $length, qlength = ". $mb->getQueryLength() . "\n";
    
    #print "$orth_protein\n\n";
    

    #print "seq2=$orth_protein\n";
    
    
    #
    # BLAST back
    # 

    # create a query file
    $fa->writeSeq($tmpfile, "ORTHOLOG", $orth_seq);

    $mb->setBlastProgram("blastn");
    $mb->setDatabaseDatabase($ARGV[0]);
    $mb->setQueryDatabase($tmpfile);

    #$mb->setVerbose(1);
    
    my $a_ref = $mb->blastallUnique;


    my $q_id = $mb->getUniqueHitName();


    if ($q_id eq $a[0]) {
	print "$a[0]\t$d_id\t$d_start\t$d_end\t$frame\n";
    }

     
    #foreach my $r (@$a_ref) {
    #print "$r->{EVALUE}\t$r->{QFROM}\t$r->{QTO}\n$r->{QSEQ}\n";
    #}

 

}


