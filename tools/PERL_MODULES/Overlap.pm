package Overlap;

use lib qw(/home/olly/PERL_MODULES);

use Sets;

require "ScanACE.pm";

# INPUT : motif, two sets of orthologous upstream seq, the list of orthologous genes
#              backgrounds for using ScanACe



sub new {
    my $self  = {};

    $self->{MOTIF_FILE}  = undef;
    
    my %h1 = ();
    $self->{BKG1}        = \%h1;

    my %h2 = ();
    $self->{BKG2}        = \%h2;

    $self->{SCORES_FILE1}        = undef;
    $self->{SCORES_FILE2}        = undef;
    $self->{PAIRS}       = [];
    
    my %h3 = ();
    $self->{INV_PAIRS}   = \%h3;

    $self->{OPT_PVALUE}    = undef;
    $self->{OPT_THRESHOLD} = undef;
    $self->{TSD_THRESHOLD} = undef;

    $self->{VERBOSE} = 0;

    $self->{NBORTH}        = undef;

    $self->{SCRIPTS_DIR} = "/home/olly/COMPARATIVE_YEAST/ORFS/BLAST_REPORTS/BINDINGMAP/SCRIPTS";

    $self->{OUTFILE} = undef;

    $self->{SHARED_ORFS} = undef;
    $self->{SHARED_POS}  = undef;
    $self->{SHARED_SEQ}  = undef;

    $self->{SPECIES1_POS} = undef;
    

    bless($self);
    return $self;
}


# returns an array of positions in the overlap set
sub getOverlapPositions {
    my ($self) = @_;

    my @a_tmp = ();
    foreach my $p1 ( @{$self->{SHARED_POS}} ) {
	foreach my $p2 ( @$p1 ) {
	    push @a_tmp, $p2;
	}
    }

    return \@a_tmp;
}


sub getSpecies1Positions {
    my ($self) = @_;

    my @a_tmp = ();
    foreach my $p1 ( @{$self->{SPECIES1_POS}} ) {
	foreach my $p2 ( @$p1 ) {
	    push @a_tmp, $p2;
	}
    }

    return \@a_tmp;
}




sub getOptimalPValue {
    
    my ($self) = @_;

    return $self->{OPT_PVALUE}; 
    
}


sub getOptimalThreshold {
    
    my ($self) = @_;
    
    return $self->{OPT_THRESHOLD}; 
    
}


#
# threshold at 2SD
#
sub get2SDThreshold {
    my ($self) = @_;
    return $self->{TSD_THRESHOLD}; 
}


#
# list of orth ORFs
#
sub setORFList {
    
    my ($self, $s_file) = @_;

    # create ORF list array
    
    open ORTH, $s_file or die "Cannot open $s_file ..\n";
    @{ $self->{PAIRS} } = <ORTH>;
    chomp @{ $self->{PAIRS} };
    close ORTH;

    $self->{NBORTH} = scalar(@{ $self->{PAIRS} });
    
    # create inverted index
    for (my $i=0; $i < $self->{NBORTH}; $i++) {
	my ($s_gene, $i_length)  = split /\t/, $self->{PAIRS}->[$i]; 
	$self->{PAIRS}->[$i]          = $s_gene;
	$self->{INV_PAIRS}->{$s_gene} = $i;

	#print "self->{INV_PAIRS}->{$s_gene} = $i;\n";

	
    }
}


# outfile containing the number of ORFs at each threshold steps
sub setOutFile {
    
    my ($self, $s_file) = @_;
    
    $self->{OUTFILE} = $s_file;
    
}

# result file to which append the 

# set the output name for the score files
sub setOutputScoreFiles {
    
    my ($self, $f1, $f2) = @_;
    
    $self->{SCORES_FILE1} = $f1;
    $self->{SCORES_FILE2} = $f2;

}


# set sequences files
sub setSequencesFiles {
    my ($self, $f1, $f2) = @_;

    $self->{SEQ_SPECIES1_FILE} = $f1;
    $self->{SEQ_SPECIES2_FILE} = $f2;
}


sub setVerbose {
    my ($self, $b) = @_;
    
    $self->{VERBOSE} = $b;
}

#
# procedure itself
#
sub calcOptimalOverlap {
    
    my ($self) = @_;

    die "Please set the names of the output score files ..\n" if (!defined( $self->{SCORES_FILE1} ) || !defined( $self->{SCORES_FILE2} ));
    die "Please provide an output file name ..\n" if (!defined($self->{OUTFILE})); 
    
    # run ScanAce for sequences 1 
    my $s1 = ScanACE->new;
    $s1->setNbMotifs(5000);
    $s1->setVerbose(1);
    $s1->setGC($self->_getGCFrequency->[0]);
    $s1->runScanACE($self->{SEQ_SPECIES1_FILE}, $self->{MOTIF_FILE});
    my $a_ref1 = $s1->getSites;
    # output results to file
    $self->_outputScoreFile($self->{SCORES_FILE1}, $a_ref1);

    # 
    $self->{TSD_THRESHOLD} = $s1->getAverage - 2 * $s1->getStdDev;

    # run ScanAce for sequences 1 
    my $s2 = ScanACE->new;
    $s2->setNbMotifs(5000);
    $s2->setGC($self->_getGCFrequency->[1]);
    $s2->runScanACE($self->{SEQ_SPECIES2_FILE}, $self->{MOTIF_FILE});
    my $a_ref2 = $s2->getSites;
    # output results to file
    $self->_outputScoreFile($self->{SCORES_FILE2}, $a_ref2);
    
    
    # LAUNCH THE COUNTS CALCULATION
    my $tmpfile   = Sets::getTempFile("/tmp/tmp");
    my $s_todo    = "$self->{SCRIPTS_DIR}/count -s1 $self->{SCORES_FILE1} -s2 $self->{SCORES_FILE2} -nborth $self->{NBORTH} -r $tmpfile >  $self->{OUTFILE}";
    print "$s_todo\n"; # if ($self->{VERBOSE} == 1);
    system(($s_todo)) == 0 or die "could not run count\n";

    #system("cat $tmpfile");

    # the previous process has create a tempo file, containing the optimal thresholds and p-values
    open IN, $tmpfile;
    my $l = <IN>; chomp $l;
    ($self->{OPT_THRESHOLD}, $self->{OPT_PVALUE}) = split /\t/, $l;
    close IN;
    
    #unlink $tmpfile;
    
}


# get both the sequence, position and orf data at the overlap 
sub buildOverlapData {

    my ($self, $t) = @_;
    #my ($i_nborfs, $f_threshold, $s_scoreFile1, $s_scoreFile2, $a_ref_pairs) = @_;

    if (!defined($t)) {
	$t = $self->{OPT_THRESHOLD};
    }
    
    
    # fill two tables of occurence
    my @a_table_pos1 = ();
    my @a_table_pos2 = ();
    
    # tables of sequences
    my @a_table_seq1 = ();
    my @a_table_seq2 = ();
    
    # initialize the tables
    for (my $i=0; $i<$self->{NBORTH}; $i++) {
	$a_table_pos1[$i] = []; 
	$a_table_pos2[$i] = []; 
	$a_table_seq1[$i] = []; 
	$a_table_seq2[$i] = []; 
    }
    
    # fill first table
    open SC1, $self->{SCORES_FILE1};
    while (my $s_line = <SC1>) {
	my ($n, $s, $p, $seq, $strand) = split /\t/, $s_line;
	if ($s >= $t) {
	    push @{$a_table_pos1[$n]}, $p;
	    push @{$a_table_seq1[$n]}, $seq;
	}	
    }
    close SC1;
    
    # fill second table
    open SC1, $self->{SCORES_FILE2};
    while (my $s_line = <SC1>) {
	my ($n, $s, $p, $seq, $strand) = split /\t/, $s_line;
	if ($s >= $t) {
	    push @{$a_table_pos2[$n]}, $p;
	    push @{$a_table_seq2[$n]}, $seq;
	}
    }
    close SC1;
    
    # output the list of shared ORF
    my @a_shared_orfs    = ();

    # list of positions in overlap
    my @a_shared_pos1    = ();
    my @a_shared_pos2    = ();

     # list of positions in species 1
    my @a_species1_pos1    = ();
    
    
    # list of sites
    my @a_shared_seq1    = ();
    my @a_shared_seq2    = ();

       
    for (my $i=0; $i<$self->{NBORTH}; $i++) {
	if ((scalar(@{$a_table_pos1[$i]}) != 0) && (scalar(@{$a_table_pos2[$i]}) != 0)) {
	    
	    push @a_shared_orfs, $self->{PAIRS}->[$i];
	    
	    # positions
	    push @a_shared_pos1, $a_table_pos1[$i];
	    push @a_shared_pos2, $a_table_pos2[$i];

	    # sites
	    push @a_shared_seq1, $a_table_seq1[$i];
	    push @a_shared_seq2, $a_table_seq2[$i];
	    
	} elsif (scalar(@{$a_table_pos1[$i]}) != 0) {

	    push @a_species1_pos1, $a_table_pos1[$i];
	    
	}
    }

    $self->{SHARED_ORFS} = \@a_shared_orfs;
    $self->{SHARED_POS} = \@a_shared_pos1;
    $self->{SHARED_SEQ} = \@a_shared_seq1;
    $self->{SPECIES1_POS} = \@a_species1_pos1;
    
    
}

#
# get the sites at the overlap
#
sub getOverlapSites {
    my ($self) = @_;
    
    my @a_tmp = ();
    foreach my $o (@{$self->{SHARED_SEQ}}) {
	foreach my $s (@{$o}) {
	    push @a_tmp, $s;
	}
    }
    return \@a_tmp;
}


#
#
#
sub getOverlapSet {
    
    my ($self) = @_;
 
    return $self->{SHARED_ORFS};

}

#
# motif that is going to be used
#
sub setMotifFile {
    my ($self, $s_file) = @_;

    $self->{MOTIF_FILE} = $s_file;
}


#
# read BKG frequencies
#
sub setBackgoundFrequencies {
    my ($self, $s_file1, $s_file2) = @_;
    
    open DB, $s_file1 or die "Could not open list of bkg frequencies $s_file1\n";
    while (my $s_line = <DB>) {
	chomp $s_line;
	my @a_tmp = split /\t/, $s_line; 
	$self->{BKG1}->{ $a_tmp[0] } = $a_tmp[1];
    }
    close DB;

    open DB, $s_file2 or die "Could not open list of bkg frequencies $s_file2\n";
    while (my $s_line = <DB>) {
	chomp $s_line;
	my @a_tmp = split /\t/, $s_line; 
	$self->{BKG2}->{ $a_tmp[0] } = $a_tmp[1];
    }
    close DB;
    
    
}


sub setBackgroundFrequenciesByHand {
    
    my ($self, $a, $c, $t, $g) = @_;
    
    $self->{BKG1}->{'C'} = $c;
    $self->{BKG1}->{'G'} = $g;    
    $self->{BKG1}->{'A'} = $a;
    $self->{BKG1}->{'T'} = $t;

    $self->{BKG2}->{'C'} = $c;
    $self->{BKG2}->{'G'} = $g;    
    $self->{BKG2}->{'A'} = $a;
    $self->{BKG2}->{'T'} = $t;

}


#
# internal method, to calc GC frequency from the BKG frequencies
#
sub _getGCFrequency {
    my ($self) = @_;
    
    return [ $self->{BKG1}->{'C'} + $self->{BKG1}->{'G'}, 
	     $self->{BKG2}->{'C'} + $self->{BKG2}->{'G'} ];
    
}


#
# output a temporary (or not) score files, easy to parse 
#  $s_scoreFile   = NAME OF THE FILE TO CREATE
#  $a_ref_motifs  = ARRAY OF MOTIFS FROM SCANACE
#  $f_threshold   = LIMITING THRESH
sub _outputScoreFile {
    my ($self, $s_scoreFile, $a_ref_motifs) = @_;
    
    open OUT, ">$s_scoreFile";
    
    foreach my $a_ref_motif (@$a_ref_motifs) {

	#print "self->{INV_PAIRS}->{$a_ref_motif->[0]}\n";
	#print "=" . $self->{INV_PAIRS}->{$a_ref_motif->[0]} . "\n";

	last if ($a_ref_motif->[4] < 0);
	print OUT $self->{INV_PAIRS}->{$a_ref_motif->[0]} . "\t";
	print OUT $a_ref_motif->[4] . "\t";
	print OUT $a_ref_motif->[1] . "\t";
	print OUT $a_ref_motif->[3] . "\t";
	print OUT $a_ref_motif->[2] . "\t";
	print OUT "\n";
    }

    close OUT;
    
}

1;

