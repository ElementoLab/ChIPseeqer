use lib qw(/home/elemento/PERL_MODULES);

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
$mb->setEvalueThreshold("1e-10");


my $file = Sets::getTempFile("toto");

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    #my @a = split /\|/, $n;

    $n =~ /(IMAGE\:\d+)/;
    
    $n = $1;

    $fa->writeSeq($file, $n, $s);
    $mb->setQueryFile($file);
    
    my $a_ref = $mb->blastallUnique;

    if (scalar(@$a_ref) == 0) {
	#print "$n\tno match\n";
	next;
    } else {
	#print scalar(@$a_ref) . " matches\n";
    }

    my $hit_name        = $mb->getUniqueHitName;
    my $hit_length      = $mb->getUniqueHitLength;
    my $query_length    = $mb->getQueryLength;

    my $a_ref_hsps      = $mb->retain_non_overlapping_blocks($a_ref, 3);

    my $a_ref_id_le     = $mb->get_indentity_and_aligned_length($a_ref_hsps);

    my $L               = $a_ref_id_le->[1]; 
    my $I               = $a_ref_id_le->[0];
    
    if (($I > 0.95) && ($L > 100)) {
	print "$n\t$hit_name\t$L\t$I\n";
    } else {
	#print "$n\tno match\n";
    }

    
}

unlink $file;
