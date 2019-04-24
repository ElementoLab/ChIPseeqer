use lib qw(/home/olly/PERL_MODULES);

use lalign;
use Fasta;
use Sets;

use File::Copy;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);
my $a_ref = $fa->nextSeq();
my ($n1, $s1) = @$a_ref;
$fa->dispose;

my $fa = Fasta->new;
$fa->setFile($ARGV[1]);


my @POS = ();
#my $i = 1;

open OUTFOOT, ">$ARGV[0].foot";

#
# go thru all the sequences, align them one by one
#

my $imax = 0;

print "Start alignments ..\n";
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
 
    my @infos = split /\|/, $n; 

    my $i     = $infos[2];     

    if ($i > $imax) {
	$imax = $i;
    }
   
    # exec lalign
    my $file = Sets::getTempFile("toto");

    $fa->writeSeq($file, $n, $s);

    print " with seq$i ($n), l=" . length($s) . " bp\n";

    my $la = lalign->new;

    
    
    $la->setMaxEvalue("1e-10");    
    $la->processUsingBlast($ARGV[0], $file);

    #$la->setMaxEvalue("0.001");
    #my $nb_asked_alns = 100;
    #$la->setNumberAlignments($nb_asked_alns);
    #$la->processUsingLalign($ARGV[0], $file);

    #
    #  re-start while the number of alns  < Emax is the asked number 
    #
    my $a_ref_res = $la->getResults;
    
    my ($rn) = $n =~ /\ (.+?)$/;
    #move($la->getLalignOutputFile(), "$ARGV[0].nn.lalign");
	
    print "  got " . scalar(@$a_ref_res) . " significant alignments\n";
    
    #
    #  analyze results
    # 
    foreach my $h (@$a_ref_res) {
	
	my $a_ref_pos = $h->{CONSERVED_POS};

	foreach my $p (@$a_ref_pos) {
	    push @{ $POS[ $p ] }, $i if !Sets::in_array($i, @{ $POS[ $p ] }); # ++;
	}
	
	print OUTFOOT "$n\t$h->{DSTART}\t$h->{DEND}\t$h->{HSEQ}\n";
	
    }

    #$i++;   # increment the index of the orth sequences (assumes that the number of genomes is always the same)
    

    unlink $file;
    Sets::unlink_blast_files($file);
    #unlink $outfile;
}


close OUTFOOT;

#my $cnt = $i-1;  # number of orth sequence

my @a = split //, $s1;

open NICE, ">$ARGV[0].nice";
my $w = 100;
for (my $i=0; $i<length($s1); $i += $w) {

    #my $p = $i+1;
    
    for (my $j=$imax; $j>=0; $j--) {
	
	for (my $k=$i; $k<$i+$w; $k++) {
	
	    #
	    #  attention here, the positions are starting at +1
	    #
	    if (defined($POS[$k+1]->[$j])) {
		print NICE $POS[$k+1]->[$j];
	    } else {
		print NICE " ";
	    }
	}
	print NICE "\n";
    }

    print NICE substr($s1, $i, $w); print NICE "\n\n";

    #print $p . "\t" . $a[$i]; print " ";
    
    #print join("", @{ $POS[$i] }); print "\n";
}
close NICE;



open OUT, ">$ARGV[0].cons";
for (my $i=0; $i<length($s1); $i += 1) {
    print OUT $a[$i]; print OUT "\t";
    
    #
    #  attention here, the positions are starting at +1
    #
    if ($POS[ $i+1 ] > 0) {
	print OUT join("", @{ $POS[$i+1] }); print OUT "\n";
    } else {
	print OUT "\n";
    }
}
close OUT;

   


    #if ($POS[ $i ] > 0) {
#	print $POS[ $i ];
#    } else {
#	print " ";
#    }
#    print "\n";
    #print "\n" if (($i % 50) == 0);
#}
#print "\n";

#for (my $i=0; $i<length($s1); $i++) {
#    print $a[$i];
    #print "\n" if (($i % 50) == 0);
#}
#print "\n";
