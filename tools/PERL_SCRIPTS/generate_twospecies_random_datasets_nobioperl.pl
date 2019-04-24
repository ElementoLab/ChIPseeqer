#
#  generate random datasets from two related species 
#  INPUT : set1, set2
#  ALGO  : foreach 2 X seq, measure divergence, learn 1M, recreate seq1, evolve to seq2
#  OUTPUT : random set1 and set2, same divergence as set1 and set2

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Sets;
use Markov;
use strict;
use Getopt::Long;
use Table;
use ClustalW;

my $s_fasta1    = undef;
my $s_fasta2    = undef;
my $s_only_orf  = undef;
my $s_after_orf = undef;
my $denovo      = 0;
my $outdir      = undef;
my $outalndir   = undef;
my $outmat      = undef;
my $freqfile    = undef;

GetOptions ('fasta1=s'     => \$s_fasta1,
            'fasta2=s'     => \$s_fasta2,
	    'only_orf=s'   => \$s_only_orf,
	    'after_orf=s'  => \$s_after_orf,
	    'denovo=s'     => \$denovo,
	    'outdir=s'     => \$outdir,
	    'outalndir=s'     => \$outalndir,
	    'outmat=s'     => \$outmat,
	    'freqfile=s'   => \$freqfile);

my $f1 = Fasta->new;
my $f2 = Fasta->new;

$f1->setFile($s_fasta1);
$f2->setFile($s_fasta2);

my $s_fasta1_out = undef;
my $s_fasta2_out = undef;

if ($outdir) {

    $s_fasta1_out = "$outdir/" . Sets::filename($s_fasta1);
    $s_fasta2_out = "$outdir/" . Sets::filename($s_fasta2);

} else {
    
    $s_fasta1_out = $s_fasta1;
    $s_fasta2_out = $s_fasta2;
}

if ($denovo == 1) {
    open OUT1, ">$s_fasta1_out.random";
}

open OUT2, ">$s_fasta2_out.random";

my $b_after_orf = 0;

my $h_ref_freqs = 0;

if (defined($freqfile)) {

    # load the entire frequency file
    my $ta = Table->new;
    $ta->loadFile($freqfile);
    $h_ref_freqs = $ta->getIndex(0);
        
}


if (defined($outmat)) {
    open OUTM, ">$outmat";
}

while (1) {
    my $a_ref1 = $f1->nextSeq;
    my $a_ref2 = $f2->nextSeq;
    
    last if (!$a_ref1 || !$a_ref2);

    my ($name1, $seq1) = @{$a_ref1};
    my ($name2, $seq2) = @{$a_ref2};

    my $len1 = length($seq1);
    my $len2 = length($seq2);


    next if ($s_only_orf && ($s_only_orf ne $name1));
    
    if ($name1 eq $s_after_orf) {
	$b_after_orf = 1;
    }
    next if ($s_after_orf && ($b_after_orf == 0));

    my %freq = ();   

    if (!defined($freqfile)) {

	#print "X" x 100;
	#print "\n";
	# global alignment (CLUSTALW)

	my $cl = ClustalW->new;
	
	$cl->setSequencesToAlign([ $name1 . "_1", $name2 . "_2" ], 
				 [ $seq1        , $seq2         ]); 
	$cl->run;
	my $a_ref_aln = $cl->getInputOrderedSeqArray;

	my $aln1 = shift @$a_ref_aln; 
	my $aln2 = shift @$a_ref_aln;

	
	# calculate the frequency of change, for each letter
	my @a1 = split //, $aln1;
	my @a2 = split //, $aln2;
	my $l  = length($aln1);
	
	my $l_actual = 0;
	for (my $i=0; $i<$l; $i++) {
	    next if (($a1[$i] =~ /[\-\.N]/) || ($a2[$i] =~ /[\-\.N]/)) ;
	    $freq{ $a1[$i] } { $a2[$i] } ++;
	    $l_actual ++;
	}
	foreach my $k1 (keys(%freq)) {
	    my $sum = 0.0;
	    foreach my $k2 (keys(%{$freq{$k1}})) {
		$freq{ $k1 } { $k2 } /= $l_actual;
		$sum += $freq{ $k1 } { $k2 };
	    }
	    
	    foreach my $k2 (keys(%{$freq{$k1}})) {
		$freq{ $k1 } { $k2 } /= $sum;
	    }
	}
	
	if (defined($outalndir)) {
	    
	    open  OUTA, ">$outalndir/$name1.seq";
	    print OUTA ">$name1\n" . $aln1 . "\n\n";
	    print OUTA ">$name2\n" . $aln2 . "\n\n";
	    close OUTA;
	    
	}


    } else {

	my @a = @{ $h_ref_freqs->{$name1} };

	shift @a;
	foreach my $r (@a) {
	    my @b = split /\//, $r;
	    $freq{ $b[0] } { $b[1] } = $b[2]; 
	    
	    #print "freq{ $b[0] } { $b[1] } = $b[2]; \n";
	}
	


    }

    #print "Subst matrix =\n";
    #foreach my $k1 (keys(%freq)) {
    #	foreach my $k2 (keys(%{$freq{$k1}})) {
    #	    print "freq{ $k1 } { $k2 } = " . $freq{ $k1 } { $k2 } . "\n";
    #	}
    #}
    

    if (defined($outmat)) {
	
	my $str = "$name1\t";
	foreach my $k1 (keys(%freq)) {
	    foreach my $k2 (keys(%{$freq{$k1}})) {
		$str .= "$k1/$k2/" . $freq{ $k1 } { $k2 } . "\t";
	    }
	}
	chop $str;

	print OUTM "$str\n";
	
	
    }


    # learn Markov model
    
    my $newseq1 = undef;

    if (($denovo == 1) && ($seq1 =~ /[ACGT]/)) {
	
	my $o_markov = Markov->new();
	$o_markov->calcFrequenciesFromSeq($seq1);
	$newseq1 = $o_markov->generate1rstOrder(length($seq1));
	
	#
	# add the Ns pattern
	#
	my @a = split //, $seq1;
	my @b = split //, $newseq1;
	for (my $i=0; $i<$len1; $i++) {
	    if ($a[$i] !~ /[ACGT]/) {
		$b[$i] = 'N';
	    }
	}
	$newseq1 = join("", @b);
	
	
    } else {
	$newseq1 = $seq1; 
    }
    
    
    
	
    # evolve the sequence
    my @a = split //, $newseq1;
    
    my $newseq2 = "";
    foreach my $n (@a) {
	    
	my $c = Sets::generateRandomSymbol($freq{$n});
	
	if ($c eq "") {
	    $c = 'N';
	}
	$newseq2 .= $c;
    }
    

    if (length($newseq1) < $len1) {
	$newseq1 = 'N' x $len1;
	
	my $len = length($newseq1);
	die "Pb for $newseq1 (should have length $len1, has $len)\n";
    }

    if (length($newseq2) < $len1) {
	#$newseq2 = 'N' x $len2;
        
        my $len = length($newseq2);
        die "Pb for $newseq2 (should have length $len2, has $len)\n";
    }
    
    if ($denovo == 1) {
	print OUT1 ">$name1\n$newseq1\n\n";
    }

    print OUT2 ">$name2\n$newseq2\n\n";
}


close OUT1;
close OUT2;
if (defined($outmat)) {
    close OUTM;
}
