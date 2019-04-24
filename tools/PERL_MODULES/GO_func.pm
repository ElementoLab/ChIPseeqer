package GO_func;
use GDBM_File;
use Hypergeom;
use Sets;

sub new {
    my $self  = {};
    
    $self->{DB_DIR} = "/home/olly/PERL_MODULES/FUNCTIONS/GO";
    
    my %t1 = ();
    my %t2 = (); 
    my %t3 = ();
    
    $self->{DB_NBFUNC_TEXT} = \%t1;
    $self->{DB_NBFUNC_NUMB} = \%t2;
    $self->{DB_ORF_NBFUNC}  = \%t3;

    #  BY DEFAULT, OPEN READ ONLY
    # takes as input a set of ORFs, and look up the most significant MIPS functions in that 
    tie(%{$self->{DB_NBFUNC_TEXT}}, 'GDBM_File', "$self->{DB_DIR}/nbfunc_text.db", &GDBM_READER, 0644);
    tie(%{$self->{DB_NBFUNC_NUMB}}, 'GDBM_File', "$self->{DB_DIR}/nbfunc_numb.db", &GDBM_READER, 0644);
    tie(%{$self->{DB_ORF_NBFUNC}},  'GDBM_File', "$self->{DB_DIR}/orf_nbfunc.db",  &GDBM_READER, 0644);
    
    #print $self->{DB_NBFUNC_NUMB};
    #print "\n";

    $self->{ORF_SET} = [];
    
    $self->{TOTALNBORFS} = 6200; #5538;  # Kellis et al number
    
    $self->{PVALUE_THRESHOLD} = undef; 

    $self->{VERBOSE} = 0;

    bless($self); 
    return $self;
    
}

sub setVerbose {

    my ($self, $i) = @_;
    $self->{VERBOSE}    = $i;
    
}


#
# GO, GO_SLIM, MIPS 
#
sub setSource {
    my ($self, $src, $spe) = @_;
    
    if (!$spe) {
	$spe = "YEAST";
    }

    if ($src eq "GO") {
	$self->{DB_DIR} = "/home/olly/PERL_MODULES/FUNCTIONS/GO/$spe";
    } elsif ($src eq "GO_SLIM") {
    	$self->{DB_DIR} = "/home/olly/PERL_MODULES/FUNCTIONS/GO_SLIM/$spe";
    } elsif ($src eq "MIPS") {
	$self->{DB_DIR} = "/home/olly/PERL_MODULES/FUNCTIONS/MIPS/$spe";
    }

    # close the ties
    untie(%{$self->{DB_NBFUNC_TEXT}});
    untie(%{$self->{DB_NBFUNC_NUMB}});
    untie(%{$self->{DB_ORF_NBFUNC}});
    
    # reopen them 
    tie(%{$self->{DB_NBFUNC_TEXT}}, 'GDBM_File', "$self->{DB_DIR}/nbfunc_text.db", &GDBM_READER, 0644);
    tie(%{$self->{DB_NBFUNC_NUMB}}, 'GDBM_File', "$self->{DB_DIR}/nbfunc_numb.db", &GDBM_READER, 0644);
    tie(%{$self->{DB_ORF_NBFUNC}},  'GDBM_File', "$self->{DB_DIR}/orf_nbfunc.db",  &GDBM_READER, 0644);

}


#
# GO, GO_SLIM, MIPS 
#
sub setDBDir {
    my ($self, $src) = @_;
    
    $self->{DB_DIR} = $src;

    # close the ties
    untie(%{$self->{DB_NBFUNC_TEXT}});
    untie(%{$self->{DB_NBFUNC_NUMB}});
    untie(%{$self->{DB_ORF_NBFUNC}});
    
    # reopen them 
    tie(%{$self->{DB_NBFUNC_TEXT}}, 'GDBM_File', "$self->{DB_DIR}/nbfunc_text.db", &GDBM_READER, 0644);
    tie(%{$self->{DB_NBFUNC_NUMB}}, 'GDBM_File', "$self->{DB_DIR}/nbfunc_numb.db", &GDBM_READER, 0644);
    tie(%{$self->{DB_ORF_NBFUNC}},  'GDBM_File', "$self->{DB_DIR}/orf_nbfunc.db",  &GDBM_READER, 0644);

}




sub readORFset {

    my ($self, $s_infile) = @_;
     
    $self->{ORF_SET} = Sets::readSet($s_infile);
    
}


#
#    fucking worm genes
#
sub renameWormGenes {
    my ($self) = @_;
    
    foreach my $o (@{$self->{ORF_SET}}) {
	
	$o =~ s/\.1$//;
	
    }

    Sets::printSet($a_ref) if ($self->{VERBOSE});


}


#
# 
# 
sub setORFset {
    my ($self, $a_ref) = @_;
    $self->{ORF_SET} = $a_ref;
    
    print "Here are the modified names :\n" if ($self->{VERBOSE});
    Sets::printSet($a_ref) if ($self->{VERBOSE});
}





#
#
sub setTotalNbORFS {
     my ($self, $n) = @_;
     
     $self->{TOTALNBORFS} = $n;    
}


#
# 
#
sub setPvalueThreshold {
    
    my ($self, $p) = @_;
     
    $self->{PVALUE_THRESHOLD} = $p;  
    
}


sub getAllAnnotatedORFs {
    
    my ($self) = @_;

    my @a = keys(%{ $self->{DB_ORF_NBFUNC} });

    return \@a;
    
    
}

#
#  return ID, texte
#
sub getAllCategories {    
    my ($self) = @_;
    return Sets::getArrayKVFromHash($self->{DB_NBFUNC_TEXT});
}


#
# get an array MIPS categories for one specific gene
#
sub getMIPScategoryIDs {
    my ($self, $orf) = @_;

    # get a set of functions
    my $s = $self->{DB_ORF_NBFUNC}->{$orf};

    # get an array of functions
    my @a = split /\|/, $s;

    return \@a;
}





#
# get an array MIPS categories for one specific gene
#
sub getMIPScategory {
    my ($self, $orf, $type) = @_;

    # get a set of functions
    my $s = $self->{DB_ORF_NBFUNC}->{$orf};

    # get an array of functions
    my @a = split /\|/, $s;

    my @a_tmp = ();

    #print keys(%{$self->{DB_NBFUNC_TEXT}});

    # count the functions
    foreach my $f (@a) {
	#print "$f\n";
	
	next if ($f == 214);
	
	#print "toto=" . $self->{DB_NBFUNC_TEXT}->{$f};
	if (defined($type)) { 
	    push @a_tmp, $f;
	} else {
	    push @a_tmp, $self->{DB_NBFUNC_TEXT}->{$f};
	}
	
    }

    return \@a_tmp;
}

#
# get the functional enrichment of an array of ORFS
#
sub getFunctionalContent {
    
    my ($self) = @_;
    
    if (scalar(@{$self->{ORF_SET}}) == 0) {
	die "Please provide a non-empty array\n";
    }
    
    my %h_cnt = ();
    my $i_cnt = 0;
    
    foreach my $o (@{$self->{ORF_SET}}) {

	my $a_ref_tmp  = $self->getMIPScategory($o, 2);
	print "$o" if ($self->{VERBOSE});
	# count the functions
	foreach my $f (@$a_ref_tmp) {
	    print "\t" . $self->{DB_NBFUNC_TEXT}->{$f} . "\n" if ($self->{VERBOSE});
	
	    $h_cnt{$f} ++;
	}
	
	print "\n"  if ($self->{VERBOSE});
	
	
	$i_cnt++;
    }
    # now traverse the non-zero categories and calc a p-value
    
    my @res = ();

    my $s2 = scalar(@{$self->{ORF_SET}});
    
    reset(%h_cnt);
    while (my ($f,$i) = each(%h_cnt)) {
    
	# get the number of gene annotated with 4f
	my $s1 = $self->{DB_NBFUNC_NUMB}->{$f};

	#print "$i,$s1,$s2,$self->{TOTALNBORFS}\n";

	my $p = Hypergeom::cumhyper($i,$s1,$s2,$self->{TOTALNBORFS});

	
	my %h_tmp = (PVALUE => $p, OVERLAP => $i, S1 => $s1, S2 => $s2, TOTAL => $s2, TEXT => $self->{DB_NBFUNC_TEXT}->{$f}, NUM => $f);

	# do nothing if p-vqlue is above threshold
	next if (defined($self->{PVALUE_THRESHOLD}) && ( $p > $self->{PVALUE_THRESHOLD} )); 
	
	push @res, \%h_tmp;
    
    }

    
    my @res_bis = sort {$a->{PVALUE} <=> $b->{PVALUE}} @res;

    return \@res_bis;
}


#
# get functional content for a specific CATEGORY
#
sub getSpecificEnrichment {
    
    my ($self, $cat) = @_;

    my $tmp;

    if (defined($self->{PVALUE_THRESHOLD})) {
	$tmp = $self->{PVALUE_THRESHOLD};
	
	$self->{PVALUE_THRESHOLD} = undef;
	
    } 

    foreach my $r (@{$self->getFunctionalContent}) {
	
	if ($r->{NUM} == $cat) {
	    $self->{PVALUE_THRESHOLD} = $tmp;
	    return $r; 
	}
    }
    
    $self->{PVALUE_THRESHOLD} = $tmp;    

    return ();
}


#
#  get all genes annotated with a given function
#
sub getAllOrfsAnnotatedWithFunction {
    
    my ($self, $f) = @_;

    my $a_ref_allorfs = $self->getAllAnnotatedORFs();

    my @a = ();

    #print "Found " . scalar(@$a_ref_allorfs) . "\n";

    foreach my $o (@$a_ref_allorfs) {
	
	my $a_ref_cat = $self->getMIPScategoryIDs($o);
	
	#print "For $o, found " . scalar(@$a_ref_cat) . " cat\n";
	
	#print join("\t", @$a_ref_cat) . "\n";
	
	#print "Is $f in there ?\n";

	if (Sets::in_array($f, @$a_ref_cat)) {
	    
	    
	    #print "in array !\n";
	    
	    my @a_tmp = ($o, 1);

	    push @a, \@a_tmp;
	} else {

	    my @a_tmp = ($o, 0);

	    push @a, \@a_tmp;
	}
    }

    return \@a;
}



#
#   open data bases for writuing
#
sub closeDatabasesForWriting {
    
    my ($self) = @_;

    # close the ties
    untie(%{$self->{DB_NBFUNC_TEXT}});
    untie(%{$self->{DB_NBFUNC_NUMB}});
    untie(%{$self->{DB_ORF_NBFUNC}});
    
    # reopen them 
    tie(%{$self->{DB_NBFUNC_TEXT}}, 'GDBM_File', "$self->{DB_DIR}/nbfunc_text.db", &GDBM_READER, 0644);
    tie(%{$self->{DB_NBFUNC_NUMB}}, 'GDBM_File', "$self->{DB_DIR}/nbfunc_numb.db", &GDBM_READER, 0644);
    tie(%{$self->{DB_ORF_NBFUNC}},  'GDBM_File', "$self->{DB_DIR}/orf_nbfunc.db",  &GDBM_READER, 0644);
}


sub resetDatabases {
    my ($self) = @_;

    # close the ties
    untie(%{$self->{DB_NBFUNC_TEXT}});
    untie(%{$self->{DB_NBFUNC_NUMB}});
    untie(%{$self->{DB_ORF_NBFUNC}});
    
    unlink "$self->{DB_DIR}/nbfunc_text.db";
    unlink "$self->{DB_DIR}/nbfunc_numb.db";
    unlink "$self->{DB_DIR}/orf_nbfunc.db";
    
}



#
#   close database for writting, open for reading
#
sub openDatabasesForWriting {
    
    my ($self) = @_;

    # close the ties
    untie(%{$self->{DB_NBFUNC_TEXT}});
    untie(%{$self->{DB_NBFUNC_NUMB}});
    untie(%{$self->{DB_ORF_NBFUNC}});
    
    # reopen them 
    tie(%{$self->{DB_NBFUNC_TEXT}}, 'GDBM_File', "$self->{DB_DIR}/nbfunc_text.db", &GDBM_WRCREAT, 0644);
    tie(%{$self->{DB_NBFUNC_NUMB}}, 'GDBM_File', "$self->{DB_DIR}/nbfunc_numb.db", &GDBM_WRCREAT, 0644);
    tie(%{$self->{DB_ORF_NBFUNC}},  'GDBM_File', "$self->{DB_DIR}/orf_nbfunc.db",  &GDBM_WRCREAT, 0644);
}


#
#  add a NBFUNC => TEXT
#
sub add_NBFUNC_TEXT {    
    my ($self, $n, $t) = @_;    
    $self->{DB_NBFUNC_TEXT}->{$n} = $t;    
}


#
#  add a NBFUNC => NUMB
#
sub add_NBFUNC_NUMB {    
    my ($self, $n, $nb) = @_;    
    $self->{DB_NBFUNC_NUMB}->{$n} += $nb;    
}


#
#  add a ORF => NBFUNC (get the current one, etc)
#
sub add_ORF_NBFUNC {    
    my ($self, $o, $n) = @_;    
    my $s = $self->{DB_ORF_NBFUNC}->{$o};
    my @a = split /\|/, $s;
    push @a, $n;
    $s = join("|", @a);
    $self->{DB_ORF_NBFUNC}->{$o} = $s;    
}


1;
