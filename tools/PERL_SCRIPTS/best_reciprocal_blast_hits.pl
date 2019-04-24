# input : set of proteins, genome
#use lib qw(/home/olly/PERL_MODULES);
BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use MyBlast;
use Fasta;
use Sets;

my $mb = MyBlast->new;
$mb->setBlastProgram("tblastn");
$mb->setDatabaseDatabase($ARGV[1]);
my $tmpfile1 = Sets::getTempFile("/tmp/blast.1");
my $tmpfile2 = Sets::getTempFile("/tmp/blast.2");

$mb->setQueryDatabase($tmpfile1);
$mb->setEvalueThreshold("1e-5");
$mb->setNbProcessors(2);
#$mb->setMismatchWeight(-);
#$mb->setVerbose(1);

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $qlen = undef;

while (my $a_ref = $fa->nextSeq) {
    
    my ($name, $seq) = @$a_ref;

    my $qlen = length($seq);

    my @a = split / /, $name;

    # create a query file
    $fa->writeSeq($tmpfile1, $name, $seq);

    # do the blast
    #$mb->setVerbose(1);
   
    $mb->setBlastProgram("tblastn");
    $mb->setDatabaseDatabase($ARGV[1]);
    $mb->setQueryDatabase($tmpfile1);

    my $a_ref = $mb->blastallUnique;
    
    next if (scalar(@$a_ref) == 0);

    # put the pieces back together, from high scoring to low scoring
    
    my @a_pieces = ();
    
    foreach my $r (@$a_ref) {

	# the current HSP is compatible if it does not overlap too much with any current piece	
	my $overlap = 0;
	foreach my $p (@a_pieces) {
	    $overlap = Sets::getSequencesOverlap($r->{QFROM}, $r->{QTO}, $p->{QFROM}, $p->{QTO});
	    last if ($overlap > 30);
	}

	if ($overlap < 30) {
	    push @a_pieces, $r;
	}
    }

    #
    # order the remanining pieces
    #
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

    my $orth_protein = "";

    my $i = 0;
    foreach my $r (@a_pieces) {
	
	#print "$r->{EVALUE}\t$r->{QFROM}\t$r->{QTO}\t$r->{DFROM}\t$r->{DTO}\t$r->{DFRAME}\n$r->{QSEQ}\n$r->{DSEQ}\n";  
	#print "\n";
	
	#$length += ($r->{QTO} - $r->{QFROM} + 1);
 
	$orth_protein .= $r->{DSEQ};
	
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
    
    

    $orth_protein =~ s/\-//g;
    my $length = length($orth_protein);

    
    #print "Total length aligned = $length, qlength = ". $mb->getQueryLength() . "\n";
    
    #print "$orth_protein\n\n";
    

    #print "seq2=$orth_protein\n";
    
    
    #
    # BLAST back
    # 

    # create a query file
    $fa->writeSeq($tmpfile2, "ORTHOLOG", $orth_protein);

    
    $mb->setBlastProgram("blastp");
    $mb->setDatabaseDatabase($ARGV[0]);
    $mb->setQueryDatabase($tmpfile2);

    #$mb->setVerbose(1);
    
    my $a_ref = $mb->blastallUnique;


    my $q_id = $mb->getUniqueHitName();

    if ($d_start > $d_end) {
      my $tt = $d_start;
      $d_start = $d_end;
      $d_end   = $tt;
    }

    if ($q_id eq $a[0]) {
	print "$a[0]\t$d_id\t$d_start\t$d_end\t$frame\n";
    }

     
    #foreach my $r (@$a_ref) {
    #print "$r->{EVALUE}\t$r->{QFROM}\t$r->{QTO}\n$r->{QSEQ}\n";
    #}

 

}

unlink $tmpfile1;
unlink $tmpfile2;
