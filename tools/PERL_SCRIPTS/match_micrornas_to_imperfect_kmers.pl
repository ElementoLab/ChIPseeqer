use lib qw(/home/olly/PERL_MODULES);

use Sets;
use Table;
use Fasta;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);

my $a_ref = $ta->getArray();

my $nb_rna_mismatches = 0;

my $i = 1;
foreach my $r (@$a_ref) {
    
    
    my $fa = Fasta->new;
    $fa->setFile($ARGV[0]);

    while (my $a_seq = $fa->nextSeq()) {
	
	my ($n, $s) = @$a_seq;

	#print "looking at $n\n";

	my $ss = uc($s);
	$ss =~ s/U/T/g;

	$ss = Sets::getComplement($ss);
	
	my $km = $r->[0];
	
	my %h = ();

	my $ov = Sets::getBestOverlap_singlestrand($ss, $km, \%h);

	$h{MATCH} =~ s/X//g;
	
	my $cov = Sets::getComplement($h{MATCH});

	# print "comparing $ss and $km, ov = $ov, $h{MATCH}\n";

	# get the sequence

	if (($ov >= (length($km)- $nb_rna_mismatches)) && ($h{RNAMATCH} <= $nb_rna_mismatches) && (length($h{MATCH}) == length($km))) {


	    # match in the transformed miRNA
	    my $lcov    =  Sets::getComplement($cov);

	    # match in the actual miRNA
	    $lcov =  lc($cov);
	    $lcov    =~ s/t/u/g;
	    
	    #print "$lcov -> $cov\n";
	    
	    $s      =~ s/$lcov/$cov/g;
	    
	    my $st = ($ov==length($km)?"":"*");

	    print "$i:$st $r->[0]\t$n\t$s\n";
	}
	
	#<STDIN>;
    }

    $i++;
}
