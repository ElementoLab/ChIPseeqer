#
#  generate random datasets from two related species 
#  INPUT : set1, set2
#  ALGO  : foreach 2 X seq, measure divergence, learn 1M, recreate seq1, evolve to seq2
#  OUTPUT : random set1 and set2, same divergence as set1 and set2


use lib qw(/home/olly/PERL_MODULES /home/olly/usr/lib/perl5/site_perl/5.6.1);
use Fasta;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::SimpleAlign;
use Sets;
use Markov;
use strict;
use Getopt::Long;


my $s_fasta1    = undef;
my $s_fasta2    = undef;
my $s_only_orf = undef;
my $s_after_orf = undef;

GetOptions ('fasta1=s'     => \$s_fasta1,
            'fasta2=s'     => \$s_fasta2,
	    'only_orf=s'  => \$s_only_orf,
	    'after_orf=s'  => \$s_after_orf);

my $factory = Bio::Tools::Run::Alignment::Clustalw->new();
my $ktuple = 4;
$factory->ktuple($ktuple);  # change the parameter before executing



my $f1 = Fasta->new;
my $f2 = Fasta->new;

$f1->setFile($s_fasta1);
$f2->setFile($s_fasta2);



open OUT1, ">$s_fasta1.sm";
open OUT2, ">$s_fasta2.sm";

my $b_after_orf = 0;

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

   

    print "X" x 100;
    print "\n";
    # global alignment (CLUSTALW)
    my $o_seq1 = Bio::Seq->new(-seq => $seq1, -id => $name1 . "_1");
    my $o_seq2 = Bio::Seq->new(-seq => $seq2, -id => $name2 . "_2");
    my @a_tmp  = ($o_seq1, $o_seq2);
    my $aln    = $factory->align(\@a_tmp);
    print $aln->percentage_identity();
    print "\n";

    # calculate the frequency of change, for each letter
    my ($o_aln1, $o_aln2)  = $aln->each_seq;    
    my @a1 = split //, $o_aln1->seq;
    my @a2 = split //, $o_aln2->seq;
    my $l  = length($o_aln1->seq);
    my %freq = ();
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
    
    print "Subst matrix =\n";
    foreach my $k1 (keys(%freq)) {
	foreach my $k2 (keys(%{$freq{$k1}})) {
	    print "freq{ $k1 } { $k2 } = " . $freq{ $k1 } { $k2 } . "\n";
	}
    }
    

    #
    #  DIFFERENCE !
    #
    my $newseq1 = $seq1;  
    
    # evolve the sequence
    my @a = split //, $newseq1;

    my $newseq2 = "";
    #print "$s_new\n\n";
    foreach my $n (@a) {
	$newseq2 .= Sets::generateRandomSymbol($freq{$n});
    }

    if (length($newseq1) < $len1) {
	$newseq1 = 'N' x $len1;
    }

    if (length($newseq2) < $len2) {
	$newseq2 = 'N' x $len2;
    }
    
    
    print OUT1 ">$name1\n$newseq1\n\n";
    print OUT2 ">$name2\n$newseq2\n\n";


    # control
    
    # global alignment (CLUSTALW)
    #my $o_seq1 = Bio::Seq->new(-seq => $newseq1, -id => $name1 . "_1");
    #my $o_seq2 = Bio::Seq->new(-seq => $newseq2, -id => $name2 . "_2");
    #my @a_tmp  = ($o_seq1, $o_seq2);
    #my $aln    = $factory->align(\@a_tmp);

    #print $aln->percentage_identity();
    #print "\n";
 
}


close OUT1;
close OUT2;
