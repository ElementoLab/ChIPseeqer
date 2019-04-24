package MyBlast;
use XML::DOM;
use Sequence;
use strict;
#use DataFiles;

#my $df = DataFiles->new;

sub new {
    my ($self) = {};
    
    $self->{BLAST_DIR}        = "/Users/olivier/PERL_MODULES/PROGRAMS/BLAST/bin"; #$df->get("BLASTDIR");
    if ($ENV{BLASTDIR} ne "") {
      $self->{BLAST_DIR} = $ENV{BLASTDIR};
    }	
    if ((! -e "$self->{BLAST_DIR}/blastall") && (! -e "$self->{BLAST_DIR}/blastall.exe")) {
      die "Cannot find blastall or blastall.exe in $self->{BLAST_DIR}. Please fix DataFiles.pm\n";
    }
    $self->{DATABASE}         = undef;
    $self->{QUERY}            = undef;
    $self->{VERBOSE}          = 0;
    $self->{EVALUE_THRESHOLD} = undef;
    $self->{STORE}            = 1;
    $self->{QUERY_LENGTH}     = undef;
    $self->{HIT_LENGTH}       = undef;
    $self->{HIT_NAME}         = undef;
    $self->{BLAST_PROGRAM}    = "blastn";
    $self->{NB_PROCESSORS}    = undef;
    $self->{CLINE}            = undef;
    $self->{EXITONCRASH}      = 0;
    $self->{CRASHED}          = 0;
    $self->{LOG}              = 0;
    $self->{RAWOUTPUT}        = undef;
    $self->{UNAME}            = `uname`; $self->{UNAME} =~ s/\n//g;
    bless($self);
    return $self;
}

sub DESTROY {
    my ($self) = @_;

    close LOG;
}


sub setExitOnCrash {
  my ($self, $e) = @_;
  $self->{EXITONCRASH} = $e;
}


sub crashed {
  my ($self) = @_;
  return $self->{CRASHED};
}


sub getRawOutput {
  my ($self) = @_;
  return $self->{RAWOUTPUT};
}

sub log {

    my ($self, $f) = @_;
    $self->{LOG} = $f;

    open LOG, ">log.txt" if ($f == 1);
}


sub setNbProcessors {
    my ($self, $f) = @_;
     $self->{NB_PROCESSORS} = $f;
}


sub setBlastProgram {
    
    my ($self, $f) = @_;
    $self->{BLAST_PROGRAM} = $f;
    
}

sub setBlastDir {
    my ($self, $f) = @_;
    $self->{BLAST_DIR} = $f;
}

sub setVerbose {
    my ($self, $f) = @_;
    $self->{VERBOSE} = $f;
}

sub setFilter {
    my ($self, $f) = @_;
    $self->{FILTER} = $f;
}

sub setMegablast {
    my ($self, $f) = @_;
    $self->{MEGABLAST} = $f;
}

sub setStoreSequences {
    my ($self, $f) = @_;
    $self->{STORE} = $f;
}


sub setDatabaseDatabase {
    my ($self, $f) = @_;
    $self->{DATABASE} = $f;
}


sub setQueryDatabase {
    my ($self, $f) = @_;
    $self->{QUERY} = $f;
}


# format then set db
sub setQueryFile {
    my ($self, $f, $t) = @_;
    
    $self->setQueryDatabase($f);
    $self->format($f, $t);
}


# format then set db
sub setDatabaseFile {
    my ($self, $f, $t) = @_;
    
    $self->setDatabaseDatabase($f);
    $self->format($f, $t);
}

#
#  add query sequence (name + sequence)
#
sub setQuerySequence {
    my ($self, $n, $s) = @_;
    
    push @{ $self->{QSEQUENCES} }, $s;
    push @{ $self->{QNAMES}     }, $n;

}


sub setEvalueThreshold {
    my ($self, $e) = @_;
    
    $self->{EVALUE_THRESHOLD} = $e;

}

sub addDatabaseSequence {
    
    
}



sub format {
    my ($self, $f, $t) = @_;
    
    my $type = "F";;
    if (defined($t)) {
	$type = $t;
    } 

    my $e = "$self->{BLAST_DIR}/formatdb -i $f -p $type -o T";
    print "$e\n" if ($self->{VERBOSE} == 1);
    
    system($e) == 0 or die "Cannot execute formatdb ..\n";
}

# blastall that only returns the top hit
# return 1 hit, with possibly multiple HSP
#  but can limit the number of hsp

sub blastallUnique {
    my ($self, $nbhsp) = @_;

    my $q_file = $self->{QUERY};
    my $d_file = $self->{DATABASE};
    
    if ($self->{UNAME} =~ /CYGWIN/) {
      $q_file = "\$(cygpath -aw \"$q_file\")";
      $d_file = "\$(cygpath -aw \"$d_file\")";
    }
    
    $self->{CRASHED} = 0;
    my $e = "$self->{BLAST_DIR}/blastall -p $self->{BLAST_PROGRAM} -i \"$q_file\" -d \"$d_file\" -v 1 -b 1 -m 7";

    
   
    $e = $self->_addOptions($e);
    

    print LOG "$e\n" if ($self->{LOG} == 1);


    print "$e\n" if ($self->{VERBOSE} == 1);
    
    $self->{CLINE} = $e;

    my $r = `$e 2> /dev/null`;  # redirect 

    #print $r;

    #if (length($r) < 50) {
    #  die "problem when executing blast with $e, got $r\n";
    #}

    return $self->_analyzeUniqueXML($r, $nbhsp);

}


#
#   blastall multiple
#
sub blastallMultiple {
    my ($self) = @_;
    
    my $q_file = $self->{QUERY};
    if ($self->{UNAME} =~ /CYGWIN/) {
      $q_file = "\$(cygpath -aw \"$q_file\")";
    }
    
    $self->{CRASHED} = 0;
    my $e = "$self->{BLAST_DIR}/blastall -p $self->{BLAST_PROGRAM} -i \"$q_file\" -d $self->{DATABASE} -m 7 ";
    
    $e = $self->_addOptions($e);
    print "$e\n" if ($self->{VERBOSE} == 1);
    $self->{CLINE} = $e;
    my $r = `$e`;

    $self->{RAWOUTPUT} = $r;

    open OUT, ">$q_file.xml" or die "Cannot open XML file for output\n";
    print OUT $r;
    close OUT;
    if (-e "$q_file.xml") {
      print STDERR "Created $q_file.xml\n";
    }


    return $self->_analyzeMultipleXML($r);

}


sub _addOptions {
    my ($self, $q) = @_;

    my $e = $q;

    if (defined($self->{EVALUE_THRESHOLD})) {
	$e .= " -e $self->{EVALUE_THRESHOLD}";
    }

    if (defined($self->{MISMATCH_WEIGHT})) {
	$e .= " -q $self->{MISMATCH_WEIGHT}";
    }
    
    if (defined($self->{GAP_OPENING})) {
	$e .= " -G $self->{GAP_OPENING}";
    }
    
    if (defined($self->{GAP_EXTENSION})) {
	$e .= " -E $self->{GAP_EXTENSION}";
    }
    
    if (defined($self->{WORD_LENGTH})) {
	$e .= " -W $self->{WORD_LENGTH}";
    }

    if (defined($self->{NB_PROCESSORS})) {
	$e .= " -a $self->{NB_PROCESSORS}";
    }

    if (defined($self->{MATCH_WEIGHT})) {
	$e .= " -r $self->{MATCH_WEIGHT}";
    }

    if (defined($self->{QUERY_STRAND})) {
      $e .= " -S $self->{QUERY_STRAND}";
    }

    if (defined($self->{FILTER}) && ($self->{FILTER} == 0)) {
      $e .= " -F F";
    }

    if (defined($self->{MEGABLAST}) && ($self->{MEGABLAST} == 1)) {
      $e .= " -n T";
    }

    return $e;
}



#
#  set the mismatch weight
#
sub setMismatchWeight {
    my ($self, $w) = @_;
    $self->{MISMATCH_WEIGHT} = $w;
}

sub setMatchWeight {
    my ($self, $w) = @_;
    $self->{MATCH_WEIGHT} = $w;
}

# idem
sub setq {
    
    my ($self, $w) = @_;
    $self->{MISMATCH_WEIGHT} = $w;
}


#
#   set gap opening weight
#
sub setGapOpening {
    my ($self, $w) = @_;
    $self->{GAP_OPENING} = $w;
    
}


#
#   set gap extension weigth
#
sub setGapExtension {
    my ($self, $w) = @_;
    $self->{GAP_EXTENSION} = $w;
    
}


#
#   set gap extension weigth
#
sub setWordLength {
    my ($self, $w) = @_;
    $self->{WORD_LENGTH} = $w;
    
}

#
#   set gap extension weigth
#
sub setQueryStrand {
    my ($self, $w) = @_;
    $self->{QUERY_STRAND} = $w;
    
}


#
# get the length of the query
#
sub getQueryLength {
    my ($self) = @_;
    return $self->{QUERY_LENGTH};
}

sub getUniqueHitLength {
    my ($self) = @_;
    return $self->{HIT_LENGTH};
}

sub getUniqueHitName {
    my ($self) = @_;
    return $self->{HIT_NAME};
}

#
#  analyze multiple hits
#
sub _analyzeMultipleXML {
    my ($self, $r, $nbhits) = @_;

    #print $r;


    # structure that contains Hits
    my @a_hits = ();

    if ($r eq "") {
      print "BLAST crashed:\n$self->{CLINE}\n"; 
      $self->{CRASHED} = 1;
      exit(0) if ($self->{EXITONCRASH} == 1);
      return \@a_hits;
    }
    

    my $parser = new XML::DOM::Parser;
    my $doc    = $parser->parse($r);
    my @a      = ();

    my $cnt    = 0;
    $self->{QUERY_LENGTH} = $doc->getElementsByTagName("BlastOutput_query-len")->item(0)->getChildNodes->item(0)->getData;
   

  
    foreach my $hit ($doc->getElementsByTagName("Hit")) {
	
	# structure that contain data related to a given hit
	my %myhit = ();

	$myhit{HIT_LENGTH}  =  $hit->getElementsByTagName("Hit_len")->item(0)->getChildNodes->item(0)->getData;
	$myhit{HIT_NAME}   =  $hit->getElementsByTagName("Hit_id")->item(0)->getChildNodes->item(0)->getData;
	$myhit{HIT_NAME}   =~ s/lcl\|//g;
	
	# array that contains all the hsps
	$myhit{HSPS} = [];

	my $Hsp_evalue_min = 100000;
	

	foreach my $hsp ($hit->getElementsByTagName("Hsp")) {
	
	    my %myhsp      = (); 
	    
	    $myhsp{EVALUE}      = $hsp->getElementsByTagName("Hsp_evalue")->item(0)->getChildNodes->item(0)->getData;
	    
	    if ($myhsp{EVALUE} < $Hsp_evalue_min) {
		$Hsp_evalue_min = $myhsp{EVALUE};
	    }
	    
	    $myhsp{QFROM}  = $hsp->getElementsByTagName("Hsp_query-from")->item(0)->getChildNodes->item(0)->getData;
	    $myhsp{QTO}    = $hsp->getElementsByTagName("Hsp_query-to")->item(0)->getChildNodes->item(0)->getData;
	    
	    if ($self->{STORE} == 1) {
		$myhsp{QSEQ}    = $hsp->getElementsByTagName("Hsp_qseq")->item(0)->getChildNodes->item(0)->getData;
	    }
	     
	    $myhsp{DFROM}    = $hsp->getElementsByTagName("Hsp_hit-from")->item(0)->getChildNodes->item(0)->getData;
	    $myhsp{DTO}      = $hsp->getElementsByTagName("Hsp_hit-to")->item(0)->getChildNodes->item(0)->getData;
	    
	    
	    $myhsp{DFRAME}   = $hsp->getElementsByTagName("Hsp_hit-frame")->item(0)->getChildNodes->item(0)->getData;
	    
	    
	    if ($self->{STORE} == 1) {
		$myhsp{DSEQ}    = $hsp->getElementsByTagName("Hsp_hseq")->item(0)->getChildNodes->item(0)->getData;
	    }

	    if ($self->{STORE} == 1) {
		$myhsp{MIDLINE}    = $hsp->getElementsByTagName("Hsp_midline")->item(0)->getChildNodes->item(0)->getData;
	    }

	    
	    $myhsp{IDENTITY}    = $hsp->getElementsByTagName("Hsp_identity")->item(0)->getChildNodes->item(0)->getData;
	    

	    $myhsp{ALIGNLEN}   = $hsp->getElementsByTagName("Hsp_align-len")->item(0)->getChildNodes->item(0)->getData;


	    push @{ $myhit{HSPS} }, \%myhsp;
	    
	    
	 }
	
	$myhit{MIN_EVALUE} = $Hsp_evalue_min;

	push @a_hits, \%myhit;
	 
     }

    
    
    $doc->dispose();
    
    return \@a_hits;

}




#
# 1 hit, multiple HSPs
#
sub _analyzeUniqueXML {
  my ($self, $r, $nbhsp) = @_;

  if ($r eq "") {
    print "BLAST crashed:\n$self->{CLINE}\n"; 
    $self->{CRASHED} = 1;
    exit(0) if ($self->{EXITONCRASH} == 1);
    return [];
  }
  

  my $parser = new XML::DOM::Parser;

  my $doc    = $parser->parse($r);

  my @a = ();

  my $cnt = 0;
  $self->{QUERY_LENGTH} = $doc->getElementsByTagName("BlastOutput_query-len")->item(0)->getChildNodes->item(0)->getData;
  my $tmp = $doc->getElementsByTagName("Hit_len")->item(0);
  if ($tmp) {
    $self->{HIT_LENGTH} = $tmp->getChildNodes->item(0)->getData;
    $self->{HIT_NAME}   = $doc->getElementsByTagName("Hit_id")->item(0)->getChildNodes->item(0)->getData;
    $self->{HIT_NAME}   =~ s/lcl\|//g;
  }

  foreach my $n ($doc->getElementsByTagName("Hsp")) {
	
	my %h_tmp = (); 
	
	$h_tmp{IDENTITY} = $n->getElementsByTagName("Hsp_identity")->item(0)->getChildNodes->item(0)->getData;
	$h_tmp{ALIGNLEN} = $n->getElementsByTagName("Hsp_align-len")->item(0)->getChildNodes->item(0)->getData;

	$h_tmp{EVALUE}   = $n->getElementsByTagName("Hsp_evalue")->item(0)->getChildNodes->item(0)->getData;
	$h_tmp{QFROM}    = $n->getElementsByTagName("Hsp_query-from")->item(0)->getChildNodes->item(0)->getData;
	$h_tmp{QTO}      = $n->getElementsByTagName("Hsp_query-to")->item(0)->getChildNodes->item(0)->getData;
	
	if ($self->{STORE} == 1) {
	    $h_tmp{QSEQ}   = $n->getElementsByTagName("Hsp_qseq")->item(0)->getChildNodes->item(0)->getData;
	}

	$h_tmp{DFROM}  = $n->getElementsByTagName("Hsp_hit-from")->item(0)->getChildNodes->item(0)->getData;
	$h_tmp{DTO}    = $n->getElementsByTagName("Hsp_hit-to")->item(0)->getChildNodes->item(0)->getData;
	
	if ($h_tmp{DFROM} > $h_tmp{DTO}) {
	  my $tmp       = $h_tmp{DFROM};
	  $h_tmp{DFROM} = $h_tmp{DTO};
	  $h_tmp{DTO}   = $tmp;
	} 

	if ($self->{BLAST_PROGRAM} ne "blastx") {
	  $h_tmp{DFRAME}  = $n->getElementsByTagName("Hsp_hit-frame")->item(0)->getChildNodes->item(0)->getData;
	}

	if ($self->{STORE} == 1) {
	    $h_tmp{DSEQ}   = $n->getElementsByTagName("Hsp_hseq")->item(0)->getChildNodes->item(0)->getData;
	}
	
	push @a, \%h_tmp;

	$cnt ++;

	last if ($cnt == $nbhsp);
	
    }
    
    $doc->dispose();

    return \@a;

}



# blastall that only returns the top hit
# return 1 hit, with possibly multiple HSP
#  but can limit the number of hsp

sub blastallMultiple_Unique {
    my ($self, $nbhsp) = @_;

    my $q_file = $self->{QUERY};
    my $d_file = $self->{DATABASE};
    
    if ($self->{UNAME} =~ /CYGWIN/) {
      $q_file = "\$(cygpath -aw \"$q_file\")";
      $d_file = "\$(cygpath -aw \"$d_file\")";
    }
    
    $self->{CRASHED} = 0;
    my $e = "$self->{BLAST_DIR}/blastall -p $self->{BLAST_PROGRAM} -i \"$q_file\" -d \"$d_file\" -v 1 -b 1 -m 7";

    
   
    $e = $self->_addOptions($e);
    

    print LOG "$e\n" if ($self->{LOG} == 1);


    print "$e\n" if ($self->{VERBOSE} == 1);
    
    $self->{CLINE} = $e;

    my $r = `$e 2> /dev/null`;  # redirect 

    open OUT, ">$q_file.out" or die "Cannot open XML file for output\n";
    print OUT $r;
    close OUT;
    if (-e "$q_file.out") {
      print STDERR "Created $q_file.out\n";
    }
    #print $r;

    #if (length($r) < 50) {
    #  die "problem when executing blast with $e, got $r\n";
    #}

    return $self->_analyzeMultiple_UniqueXML($r, $nbhsp);

}



#
# multiple runs, each with 1 hit, multiple HSPs
#
sub _analyzeMultiple_UniqueXML {
  my ($self, $r, $nbhsp) = @_;

  if ($r eq "") {
    print "BLAST crashed:\n$self->{CLINE}\n"; 
    $self->{CRASHED} = 1;
    exit(0) if ($self->{EXITONCRASH} == 1);
    return [];
  }
  
  my @a_runs = ();  # 1 run = 1 read

  my $parser = new XML::DOM::Parser;
  my $doc    = $parser->parse($r);
  
  my $cnt = 0;

  foreach my $iter ($doc->getElementsByTagName("Iteration")) {
    
    my %query = ();    
    # get length
    $query{LENGTH} = $iter->getElementsByTagName("Iteration_query-len")->item(0)->getChildNodes->item(0)->getData;
    # get ID
    $query{NAME} = $iter->getElementsByTagName("Iteration_query-def")->item(0)->getChildNodes->item(0)->getData;   
    #print "$query{NAME}";

    # get first and unique hit
    my %besthit = ();
    my $hit = $iter->getElementsByTagName("Hit")->item(0);
    if ($hit) {
      $besthit{NAME} = $hit->getElementsByTagName("Hit_id")->item(0)->getChildNodes->item(0)->getData;
      $besthit{NAME} =~ s/lcl\|//g;
      $besthit{LENGTH} = $hit->getElementsByTagName("Hit_len")->item(0)->getChildNodes->item(0)->getData;

      #print "\t$besthit{NAME}";
    } else {
      #print "\tNo match";
    }
    #print "\n";
    
    #hsps
    my @hsps = ();
    foreach my $n ($hit->getElementsByTagName("Hsp")) {
	
      my %h_tmp = (); 
      
      $h_tmp{IDENTITY} = $n->getElementsByTagName("Hsp_identity")->item(0)->getChildNodes->item(0)->getData;
      $h_tmp{ALIGNLEN} = $n->getElementsByTagName("Hsp_align-len")->item(0)->getChildNodes->item(0)->getData;
      
      $h_tmp{EVALUE}   = $n->getElementsByTagName("Hsp_evalue")->item(0)->getChildNodes->item(0)->getData;
      $h_tmp{QFROM}    = $n->getElementsByTagName("Hsp_query-from")->item(0)->getChildNodes->item(0)->getData;
      $h_tmp{QTO}      = $n->getElementsByTagName("Hsp_query-to")->item(0)->getChildNodes->item(0)->getData;
      
      $h_tmp{QSEQ}   = $n->getElementsByTagName("Hsp_qseq")->item(0)->getChildNodes->item(0)->getData;
      
      $h_tmp{DFROM}  = $n->getElementsByTagName("Hsp_hit-from")->item(0)->getChildNodes->item(0)->getData;
      $h_tmp{DTO}    = $n->getElementsByTagName("Hsp_hit-to")->item(0)->getChildNodes->item(0)->getData;
      
      if ($h_tmp{DFROM} > $h_tmp{DTO}) {
	my $tmp       = $h_tmp{DFROM};
	$h_tmp{DFROM} = $h_tmp{DTO};
	$h_tmp{DTO}   = $tmp;
      } 
      
      $h_tmp{DFRAME}  = $n->getElementsByTagName("Hsp_hit-frame")->item(0)->getChildNodes->item(0)->getData;
      $h_tmp{QFRAME}  = $n->getElementsByTagName("Hsp_query-frame")->item(0)->getChildNodes->item(0)->getData;

      $h_tmp{DSEQ}   = $n->getElementsByTagName("Hsp_hseq")->item(0)->getChildNodes->item(0)->getData;
	
      push @hsps, \%h_tmp;
    } # loop over hsps
    
    my %run = ("QUERY"=>\%query, "HIT"=>\%besthit, "HSPS"=>\@hsps);
    push @a_runs, \%run;

  } # loop over iter
  $doc->dispose();
  
  return \@a_runs;

}






#
#  the first HSP defines the anchor; the query is projected onto the best HSP, defining approximate boundaries; 
#     all HSPs that overlap the projection
#

sub getBoundariesOfQueryMatch {
    my ($self, $a_ref) = @_;

    # from the first hit, build the approx boundaries by simple projection
    #$approx_from = 
    
    
}

#
#  shortcut
#
sub getSequence {
    my ($self, $db, $name, $start, $end) = @_;
    my $s = Sequence->new;
    $s->setBlastPath( $self->{BLAST_DIR} );
    $s->setBlastDB( $db );

    if ($self->{VERBOSE} == 1) {
	$s->setVerbose(1);
    }

    my $seq = $s->getSequenceFromBlastDB($name, $start, $end);
    return $seq;
}


sub retain_non_overlapping_blocks {

    my ($self, $a_ref_hsps, $ov_max) = @_;
    
    my @a_pieces = ();
    
    foreach my $r (@$a_ref_hsps) {
    	
	my $overlap = 0;
	foreach my $p (@a_pieces) {
	    $overlap = Sets::getSequencesOverlap($r->{QFROM}, $r->{QTO}, $p->{QFROM}, $p->{QTO});
	    last if ($overlap > $ov_max);
	}
	
	if ($overlap <= $ov_max) {
	    push @a_pieces, $r;
	}

    }

    
    @a_pieces = sort { $a->{QFROM} <=> $b->{QFROM} } @a_pieces;
    
    return \@a_pieces;
    
}


sub get_query_total_aligned_length {
    
    my ($self, $a_ref_hsps) = @_;
    
    my $l1 = 0;
    foreach my $h (@$a_ref_hsps) {
	my $s1 =  $h->{"QSEQ"};             
	$s1 =~ s/\-//g;
	$l1 += length($s1);
    }

    return $l1;
}


sub getHSPdata {
  my ($self, $qseq, $dseq) = @_;

  # raw HSP len
  my $ug_qseq = $qseq;
  $ug_qseq =~ s/\-//g;
  my $ug_qseq_len = length($ug_qseq);
  
  # num matches
  my @a           = split //, $qseq;
  my @b           = split //, $dseq;
  my $len         = scalar(@a);
  my $cnt_matches = 0;
  for (my $j=0; $j<$len; $j++) {
    $cnt_matches ++ if ($a[$j] eq $b[$j]);
  }
    
  # num insertions
  my $numdel = 0;
  while ($qseq =~ /\-+/g) {
    $numdel++;
  }

  # num del
  my $numins = 0;
  while ($dseq =~ /\-+/g) {
    $numins++;
  }

  return [ $ug_qseq_len, $cnt_matches, $numins, $numdel ];
}


sub get_indentity_and_aligned_length {
    
    my ($self, $a_ref_hsps) = @_;

    my $l1 = 0;
    my $l2 = 0;
    my $gapped_seq1 = "";
    my $gapped_seq2 = "";
    
    foreach my $h (@$a_ref_hsps) {
	my $s1 =  $h->{"QSEQ"};             
	$gapped_seq1 .= $s1;
	$s1 =~ s/\-//g;
	$l1 += length($s1);
	my $s2 =  $h->{"DSEQ"};  
	$gapped_seq2 .= $s2;
	$s2 =~ s/\-//g;
	$l2 += length($s2);
    }
    
    #  compute the identity between the natching segments
    my @a           = split //, $gapped_seq1;
    my @b           = split //, $gapped_seq2;
    my $len         = scalar(@a);
    my $cnt_matches = 0;
    my $cnt_gaps1   = 0;
    my $cnt_gaps2   = 0;
    for (my $j=0; $j<$len; $j++) {
	$cnt_matches ++ if ($a[$j] eq $b[$j]);
	$cnt_gaps1   ++ if ($a[$j] eq '-');
	$cnt_gaps2   ++ if ($b[$j] eq '-');
    }
    
    #
    #   calculate the real fraction of identical residues
    #
    
    my $frac_identical = ($cnt_matches / ( length($gapped_seq1) - $cnt_gaps1 - $cnt_gaps2 ));

    
    my $L = length($gapped_seq1) - ( $cnt_gaps1 + $cnt_gaps2 ) / 2;

    return [ $frac_identical, $L ];
    
}

1;
