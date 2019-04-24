#!/usr/bin/perl

use MyBlast;
use Fasta;

my $fa = Fasta->new;

$fa->setFile($ARGV[0]);


my $mb = MyBlast->new;

$mb->setTargetDatabase("/home/olly/GENOMES/HUMAN_MOUSE/alu.n");

while (my $a_ref = $fa->nextSeq) {

    my ($name, $seq) = @{$a_ref};
	

    #print length($seq);
    
    $fa->writeSeq("tmp.fasta", $name, $seq);
    
    $mb->setQueryFile("tmp.fasta");
    
    my $a_ref = $mb->blastallUnique;
    
    foreach my $r (@$a_ref) {    
	#print "$r->{EVALUE}\t$r->{QFROM}\t$r->{QTO}\n$r->{QSEQ}\n";  

	next if ($r->{EVALUE} > 1e-20);

	if ($r->{QFROM} > $r->{QTO}) {
	    my $tmp = $r->{QTO};
	    $r->{QTO} = $r->{QFROM};
	    $r->{QFROM} = $tmp;
	}
	
	substr($seq, $r->{QFROM} - 1, $r->{QTO} - $r->{QFROM} + 1)  = 'N' x ( $r->{QTO} - $r->{QFROM} + 1);
	
	#print "$seq\n";
    }
    
    
    print ">$name\n$seq\n\n";
}
