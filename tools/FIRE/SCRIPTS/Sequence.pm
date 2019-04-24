package Sequence;
#use Bio::SearchIO;


sub new {
    my ($self)  = {};
     
    $self->{FASTA}     = undef;
    $self->{BLASTPATH} = "/home/elemento/PERL_MODULES/PROGRAMS/BLAST";

    $self->{BLASTDB}   = undef;
    $self->{NOTFOUND}  = [];
    
    $self->{NAME}      = undef;
    $self->{SEQUENCE}  = undef;
    
    # alignment variables
    $self->{SEQ1_ALIGNED} = undef;
    $self->{SEQ2_ALIGNED} = undef;
    
    # identity
    $self->{DISTANCE}    = undef;
    $self->{DISTANCE_JC} = undef;
    
    $self->{UNAME}            = `uname`; $self->{UNAME} =~ s/\n//g;
    
    $self->{USEGI}       = 0;

    # identity corrected

    $self->{VERBOSE}   = 0;
    bless($self);           # but see below
    return $self;
}

sub setFastaProgram {
    my ($self, $f) = @_;

    $self->{FASTA}     = $f;

}


sub getSequencesNotFound {
    my ($self) = @_;
    return $self->{NOTFOUND};
}


sub seq {
     my ($self, $s) = @_;
     $self->{SEQUENCE} = $s;
}


sub getSeq {
     my ($self) = @_;
     return $self->{SEQUENCE};
}

#
# added for IG annot
#
sub setSeq {
     my ($self, $n, $s) = @_;
     $self->{SEQUENCE} = $s;
     $self->{NAME}     = $n;
 }

#
# mask the seq
#
sub mask {
    my ($self, $p1, $p2) = @_;

    print "Masking from " . ($p1 - 1) . ", length " . ($p2 - $p1 + 1) . "\n" if ($self->{VERBOSE} == 1);

    substr($self->{SEQUENCE}, $p1 - 1, $p2 - $p1 + 1)  = 'N' x ( $p2 - $p1 + 1);
}


#
#  write the sequence to disk
#
sub writeToFasta {
    my ($self, $f) = @_;
	 
    open OUTF, ">$f" or die "Cannot open $f\n";
    print OUTF ">$self->{NAME}\n$self->{SEQUENCE}\n\n";
    close OUTF;
    
}

sub getComplement {
    
    my ($self, $str) = @_;
    
    my $l = length($str);
    
    my @s = split //, $str;

    my $c = "";
    for (my $i=$l-1; $i>=0; $i--) {
	
	my $d = "";
	
	if ($s[$i] eq 'A') {
	    $d = 'T';
	} elsif ($s[$i] eq 'T') {
	    $d = 'A';
	} elsif ($s[$i] eq 'G') {
	    $d = 'C';
	} elsif ($s[$i] eq 'C') {
	    $d = 'G';
	} else {
	    $d = 'N';
	}

	$c .= $d;
    }
    
    return $c;
   
}



#
#  
#
sub setVerbose {
    my ($self, $n) = @_;
    $self->{VERBOSE} = $n;
}

sub useGI {
  my ($self, $n) = @_;
  $self->{USEGI} = $n;
}

sub setBlastPath {
    
    my ($self, $s_path) = @_;

    die "Blast path incorrect\n" if (! -e "$s_path/fastacmd");  
    
    $self->{BLASTPATH} = $s_path;
       
}

sub setBlastDB {
    
    my ($self, $s_db) = @_;

    if (! -e $s_db) {
	die "$s_db does not exist ..\n";
    }

    $self->{BLASTDB} = $s_db;
       
}

sub getSequenceFromBlastDB {
    my ($self, $s_seq, $i_start, $i_end) = @_;

    die "Pb : BLAST path not defined ..\n" if (!$self->{BLASTPATH});
    die "Pb : BLAST db not defined ..\n" if (!$self->{BLASTDB});
    
    
    my $d_file = $self->{BLASTDB};
    if ($self->{UNAME} =~ /CYGWIN/) {
      $d_file = "\$(cygpath -aw \"$d_file\")";
    }

    my $seqid = undef;
    if ($self->{USEGI} == 1) {
      $seqid = "gi|$s_seq";
    } else {
      $seqid = "lcl|$s_seq";
    }	

    my $s_todo = "$self->{BLASTPATH}/fastacmd -d $d_file -s '$seqid' -L $i_start,$i_end 2> /dev/null";

    print "$s_todo\n" if ($self->{VERBOSE} == 1);

    my $s_sequence   = `$s_todo`;
    $s_sequence =~ s/^\>(.+)\n?//;
    $s_sequence =~ s/\n//g;
    $s_sequence =~ s/\r//g;


    return $s_sequence;
}

sub addSequenceToFile {
    
    my ($self, $seq, $name, $file) = @_;
    
    open OUT, ">>$file";
    
    print OUT ">$name\n$seq\n\n";

    close OUT;

}


sub writeSequencesToFile {
    
    my ($self, $a_ref, $s_file) = @_;
	
    open OUT, ">$s_file" or die "Could not create file $s_file\n";
    foreach my $o (@$a_ref) {
	my $s1 = $self->getSequenceFromBlastDB($o, 0, 0);
	
	if ($s1) {
	    print OUT ">$o\n$s1\n\n";
	} else {
	    push @{$self->{NOTFOUND}}, $o;
	}
    }
    close OUT;

}


# returns the ORFs of a given DNA sequence
sub getORFs {
    
    my ($self) = @_;
    
    # get  the positions of the ATG
    my @a_atg_pos = $self->_match_positions("ATG", $self->{SEQUENCE});

    # get the positions of the STOPs
    my @a_stop_pos = $self->_match_positions("(TA[AG]|TGA)", $self->{SEQUENCE});
    
   
    #print join("\n", @a_atg_pos) . "\n\n";
    
    #print join("\n", @a_stop_pos) . "\n";


    # array of ORF sequences
    my @a_orfs     = ();

    # foreach ATG, get the position of the first inframe STOP
    
    foreach $atg (@a_atg_pos) {

	my $good_atop = undef;

	foreach $stop (@a_stop_pos) {
	    
	    next if ($stop <= $atg);

	    $good_stop = $stop, last if ((($stop - $atg) % 3) == 0);
	    
	    
	}
	
	if ($good_stop) {
	    my %h = (S => $atg, E => $good_stop, SEQ => substr($self->{SEQUENCE}, $atg, $good_stop - $atg));
	    push @a_orfs, \%h;
	}

    }
    
    
    return \@a_orfs;

}

sub _match_positions{

    my ($self, $regexp, $sequence) = @_;
    
    

    my @positions =();
    
    while ($sequence =~/$regexp/ig){
	push (@positions, pos($sequence) - length($&) );
    }

    

    return @positions;
}


# suppose we have the start and end position of the real seq
#    we can determine easily the overlap between tzo seq
# s1 s2  e1 e2
sub getSimpleOverlap {
    
    my ($self, $s1, $e1, $s2, $e2) = @_;

    print "overlap between $s1, $e1, $s2, $e2\n";
    
    return 0 if ($s2 > $e1);
    return 0 if ($s1 > $e2);

    my @a_tmp = ($s1, $e1, $s2, $e2);
    
    @a_tmp_sorted    = sort { $a <=> $b } @a_tmp;

    #print "sorted = " . join("-", @a_tmp_sorted) . "\n";

    return $a_tmp_sorted[2] - $a_tmp_sorted[1];
    
    
}

# run a simple FASTA between two DNA sequence
sub runFasta {

    my ($self, $s_seq1, $s_seq2) = @_;

    if (!$self->{FASTA} || ! -e $self->{FASTA}) {
	die 
    }
    
    my $s_tmpstore1 = "/tmp/fasta1.tmp";
    my $s_tmpstore2 = "/tmp/fasta2.tmp";
    my $s_tmpreport = "/tmp/fastar.tmp";
    
    # create two  temporary file
    open  SEQ1,">$s_tmpstore1";
    print SEQ1 ">" . "SEQ1" . "\n" . $s_seq1 . "\n\n";
    close SEQ1;
    open  SEQ2,">$s_tmpstore2";
    print SEQ2 ">" . "SEQ2" . "\n" . $s_seq2 . "\n\n";
    close SEQ2;

    # use FASTA to align the current sequence against the reference sequence
    $s_output = `$self->{FASTA} -A -Q -b1 $s_tmpstore1 $s_tmpstore2`;

    # output the report
    open  OUT, ">$s_tmpreport";
    print OUT  $s_output;
    close OUT;

    # create a BioPerl object to handle the report
    my $searchio = new Bio::SearchIO(-format => 'fasta',
                                     -file   => $s_tmpreport);
    my $result = $searchio->next_result();
    my $hit    = $result->next_hit;
    my $hsp    = $hit->next_hsp;
    $self->{DISTANCE}    = 1 - $hsp->frac_identical;    
    $self->{DISTANCE_JC} = ($self->{DISTANCE} < 0.75 ) ? -3.0/4.0*log(1-(4.0*$self->{DISTANCE}/3.0)) : "NA";
    
    

    # get the query string
    $s_seq1 = $hsp->query_string;

    # remove the gaps
    #$s_seq1 =~ s/\-//g;
    
    $s_seq2 = $hsp->hit_string;

    #print "$s_seq1\n";
    #print $hsp->homology_string;
    #print "\n";
    #print "$s_seq2\n";
    
    #print $hsp->length('total') . "\n";


    
    
}

#
# 
#
sub getDistance {
    my ($self) = @_;
    return $self->{DISTANCE};
}


#
# 
#
sub getDistanceJC {
    my ($self) = @_;
    return $self->{DISTANCE_JC};
}



# translate a sequence

sub translate {
    my ($self) = @_;

    #die "Sequence.pm: Length of seq to translate not a multiple of 3\n" if ((length($self->{SEQUENCE}) % 3) != 0);
    
    my %CODON_TABLE = (
   TCA => 'S',TCG => 'S',TCC => 'S',TCT => 'S',
   TTT => 'F',TTC => 'F',TTA => 'L',TTG => 'L',
   TAT => 'Y',TAC => 'Y',TAA => '*',TAG => '*',
   TGT => 'C',TGC => 'C',TGA => '*',TGG => 'W',
   CTA => 'L',CTG => 'L',CTC => 'L',CTT => 'L',
   CCA => 'P',CCG => 'P',CCC => 'P',CCT => 'P',
   CAT => 'H',CAC => 'H',CAA => 'Q',CAG => 'Q',
   CGA => 'R',CGG => 'R',CGC => 'R',CGT => 'R',
   ATT => 'I',ATC => 'I',ATA => 'I',ATG => 'M',
   ACA => 'T',ACG => 'T',ACC => 'T',ACT => 'T',
   AAT => 'N',AAC => 'N',AAA => 'K',AAG => 'K',
   AGT => 'S',AGC => 'S',AGA => 'R',AGG => 'R',
   GTA => 'V',GTG => 'V',GTC => 'V',GTT => 'V',
   GCA => 'A',GCG => 'A',GCC => 'A',GCT => 'A',
   GAT => 'D',GAC => 'D',GAA => 'E',GAG => 'E',
   GGA => 'G',GGG => 'G',GGC => 'G',GGT => 'G');
    
    my @a = split //, $self->{SEQUENCE};
    my $l = length($self->{SEQUENCE});
    
    my $t = "";
    for (my $i=0; $i<$l; $i+=3) {
	my $c = $CODON_TABLE { $a[$i] . $a[$i+1]  . $a[$i+2] };
	$t .= ($c?$c:"?");
    }

    return $t;
}


# shuffle a sequence

1;

