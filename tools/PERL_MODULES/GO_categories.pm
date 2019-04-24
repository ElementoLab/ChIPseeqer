package GO_categories;
use lib qw(/home/olly/PERL_MODULES);
use Database;
use Hypergeom;
use Sets;

sub new {
    my $self  = {};
    

    $self->{DB} = Database->new;
    
    #
    # by default, connect to the human database
    #
    #$self->{DB}->connect("HUMAN");


    $self->{USER}             = undef;
    $self->{PASS}             = undef;
    $self->{HOST}             = undef;

    $self->{ORF_SET}          = [];    
    $self->{TOTALNBORFS}      = 30000; #6200; #5538;  # Kellis et al number    
    $self->{PVALUE_THRESHOLD} = undef; 
    $self->{VERBOSE}          = 0;
    $self->{ORIGINAL}         = 0;
    #
    $self->{SOURCE}           = "mysql";

    $self->{ANNOTATION_TABLE} = "GO_FULL_ANNOTATION";
    
    $self->{NBCATEGORIES}     = 0;
    
    $self->{BONFERRONI}       = 0;


    $self->{TYPE}       = undef;




    bless($self); 
    return $self;
    
}

sub setAnnotationType {

    my ($self, $t) = @_;
    $self->{TYPE} = $t;
}

sub setSource {
    my ($self, $i) = @_;
    $self->{SOURCE} = $i;

    $self->{DB}->connect("HUMAN");

}

sub setBonferroniCorrection {
    my ($self, $i) = @_;
     $self->{BONFERRONI} = $i;
     
}

#
#  
#
sub setVerbose {
    my ($self, $i) = @_;
    $self->{VERBOSE}    = $i;

    #if ($self->{VERBOSE} == 1) {
    #$self->{DB}->setVerbose($i);
    #}
}

#
#  set the ids for the connection
#
sub setID {
    my ($self, $u, $p, $h) = @_;
    
    $self->{USER} = $u;
    $self->{PASS} = $p;
    $self->{HOST} = $h;
    
}



#
#  select the species
#
sub setSpecies {
    my ($self, $s) = @_;

    if ($s ne $self->{SPECIES}) {
	$self->{SPECIES}    = $s;
	#$self->{DB}->disconnect;

	#
	#  set the id if needed
	#
	if (defined($self->{USER}) && defined($self->{PASS}) && defined($self->{HOST})) {
	    $self->{DB}->setID($self->{USER}, $self->{PASS}, $self->{HOST});
	}

	$self->{DB}->connect($s);

	if ($self->{VERBOSE} == 1) {
	    print "Connecting to database $s\n";
	}
    } 
}

sub setAnnotationTable {

    my ($self, $table) = @_;     
    $self->{ANNOTATION_TABLE} = $table;

}


#
#  
#
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

    Sets::printSet($self->{ORF_SET}) if ($self->{VERBOSE});
}


#
# 
# 
sub setORFset {
    my ($self, $a_ref) = @_;
    $self->{ORF_SET} = $a_ref;
    
    #print "Here are the modified names :\n" if ($self->{VERBOSE});
    Sets::printSet($a_ref) if ($self->{VERBOSE});
}


#
# set ENSEMBL ORF set (transform ENSEMBL -> SWISSPROT)
#
sub setEnsemblORFset {
    my ($self, $a_ref) = @_;
    
    #
    # build SQL query to get the ENSEMBL ids
    #
    my $set = Sets::SetToSQLSet($a_ref, "'");
    my $sql = "SELECT distinct SWISSPROT_ID from ENSEMBL_GENES where ENSEMBL_ID in $set and SWISSPROT_ID != '' group by ENSEMBL_ID";
    
    #print "Sql=$sql\n";
    my $a_ref_swiss = $self->{DB}->queryAllRecordsRef($sql);
    $self->{ORF_SET} = Sets::SQLRefToSet($a_ref_swiss, "SWISSPROT_ID");
    
    
    Sets::writeSetToFile($self->{ORF_SET}, "toto.txt"); 
}


sub getORFset {
    my ($self) = @_;
    
    return $self->{ORF_SET};
}


#
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


sub setOriginal {
    my ($self, $ori) = @_;
    
    $self->{ORIGINAL} = $ori;
}

#
# get an array MIPS categories for one specific gene
#
sub getGOCategories {
    my ($self, $orf) = @_;

    #
    #  prepare the mysql database
    #

    

    if ($self->{SOURCE} eq "mysql") {
	
	my $sql   = undef;
	if ($self->{SPECIES} eq "HUMAN") {
	    $sql = "SELECT distinct GOID FROM $self->{ANNOTATION_TABLE} WHERE US  = '$orf'";
	} elsif (($self->{SPECIES} eq "DROSOPHILA") || ($self->{SPECIES} eq "CIONA")) {
	    $sql = "SELECT distinct GOID FROM $self->{ANNOTATION_TABLE} WHERE UID = '$orf'";
	} elsif ($self->{SPECIES} eq "WORMS") {
	    $sql = "SELECT distinct GOID FROM $self->{ANNOTATION_TABLE} WHERE UID = '$orf'";
	}
	
	if ( $self->{ORIGINAL} == 1) {
	    $sql .= " AND ORIGINAL = 1";
	}

	if ( defined($self->{TYPE}) ) {
	    $sql .= " AND ASPECT = '$self->{TYPE}'";
	}

	if ($self->{VERBOSE} == 1) {
	    print "$sql\n";
	}
	
	my $a_ref = $self->{DB}->queryAllRecordsRef($sql);

	my @tmp = ();
	foreach my $r (@$a_ref) {
	    push @tmp, $r->{GOID};
	}
	
	
	return \@tmp;
    
    
    }  else {
	
	return $self->{DB_ORF_NBFUNC}->{$orf};
	
    }
 


    
}

#
#  special fonction for reading arbitraty ID => CAT1,CAT2 annotation
#
sub readAnnotation {
    
    my ($self, $file) = @_;

    my $ta = Table->new;
    $ta->setFile($file);
    
    while (my $r = $ta->nextRow()) {
	
	if ($r->[1] ne "") {
	    
	    my @a = split /\,/, $r->[1];

	    my $a_ref_unique = Sets::removeDuplicates(\@a);

	    $self->{DB_ORF_NBFUNC}->{$r->[0]} = $a_ref_unique;

	    foreach my $g (@$a_ref_unique) {
		$self->{DB_NBFUNC_NUMB}->{$g} ++;
	    }
	    
	}
	
    }

    $self->{NBCATEGORIES} = scalar(keys(%{ $self->{DB_NBFUNC_NUMB} }));

}



#
#  get all info related to a given GO category 
#
sub getGOSingleCategoryNumbers {
    
    my ($self, $goid) = @_;
    
    if ($self->{SOURCE} eq "mysql") {
	
	my $sql = undef;
	if ($self->{SPECIES} eq "HUMAN") {
	    $sql = "SELECT count(distinct(US))  as COUNT FROM GO_FULL_ANNOTATION where GOID = '$goid' and US like '%HUMAN'";
	} elsif (($self->{SPECIES} eq "DROSOPHILA") || ($self->{SPECIES} eq "CIONA")) {
	    $sql = "SELECT count(distinct(UID)) as COUNT FROM GO_FULL_ANNOTATION where GOID = '$goid'";
	} elsif ($self->{SPECIES} eq "WORMS") {
	    $sql = "SELECT count(distinct(UID)) as COUNT FROM GO_FULL_ANNOTATION where GOID = '$goid'";
	}

	if ( defined($self->{TYPE}) ) {
	    $sql .= " AND ASPECT = '$self->{TYPE}'";
	}

	if ($self->{VERBOSE} == 1) {
	    print "$sql\n";
	}

	my $h_ref = $self->{DB}->queryOneRecord($sql);
	
	return $h_ref->{COUNT};
    
    
    } else {
	
	return $self->{DB_NBFUNC_NUMB}->{$goid};

    }
}


sub setMaxCategory {
    my ($self, $m) = @_;
    $self->{MAXCATEGORY} = $m;
}

sub setMinCategory {
    my ($self, $m) = @_;
    $self->{MINCATEGORY} = $m;
}


#
#  get GOID, NAME for all categories
#
sub getGOCategoriesNames {
    my ($self) = @_;
    my $sql = "SELECT GOID, NAME FROM GO_CATEGORIES";

    if ($self->{VERBOSE} == 1) {
	print "$sql\n";
    }

    my $a_ref = $self->{DB}->queryAllRecordsRef($sql);    
    my %h = ();
    foreach my $r (@$a_ref) {
	if ($self->{VERBOSE} == 1) {
	    #print "$r->{GOID} => $r->{GOID} / $r->{NAME} \n";
	}
	$h{ $r->{GOID} } = $r;
    }
    return \%h;    
}

sub getAllGOcategories {
    my ($self) = @_;
    my $sql = "SELECT GOID, NAME FROM GO_CATEGORIES";

    

    my $a_ref = $self->{DB}->queryAllRecordsRef($sql);    
    my @h = ();
    foreach my $r (@$a_ref) {
	my @a = ($r->{GOID}, $r->{NAME});
	push @h, \@a;
    }
    return \@h;    
}


sub getAllGenesInCategory {
    my ($self, $goid) = @_;
    
    $sql = "SELECT distinct(UID) FROM GO_FULL_ANNOTATION where GOID = '$goid'";

    if ( defined($self->{TYPE}) ) {
	$sql .= " AND ASPECT = '$self->{TYPE}'";
    }

    my $a_ref = $self->{DB}->queryAllRecordsRef($sql);    
    my @set = ();
    foreach my $r (@$a_ref) {
	push @set, $r->{UID};
    }
    return \@set;   
}


#   
#  get the functional enrichment of an array of ORFS
#
sub getFunctionalContent {
    
    my ($self) = @_;
    
    if (scalar(@{$self->{ORF_SET}}) == 0) {
	die "Please provide a non-empty array\n";
    }
    
    my %h_cnt = ();
    my $i_cnt = 0;


    #  
    #   get GOID, NAME for all GO categories, in a hash GOID => [GOID, NAME]
    #
    my $h_ref_names = $self->getGOCategoriesNames();
    
    #
    #   for each ORF in the set, get the function with which it is annotated 
    #
    foreach my $o (@{$self->{ORF_SET}}) {
	print "Look at Gene $o\n" if ($self->{VERBOSE});
	#
	# get all the categories
	#
	my $a_ref_gofunc  = $self->getGOCategories($o);

	#
	# count the functions
	#
	foreach my $f (@$a_ref_gofunc) {
	    print " func $f $h_ref_names->{$f}->{NAME}\n" if ($self->{VERBOSE});
	    $h_cnt{$f} ++;
	}
	
	$i_cnt++;
    }

    $self->{NBCATEGORIES} = scalar(keys(%h_cnt));

    #
    # now traverse the non-zero categories and calc a p-value
    #
    my @a_res = ();

    my $s2 = scalar(@{$self->{ORF_SET}});

    reset(%h_cnt);
    while (my ($f,$i) = each(%h_cnt)) {
    	
	#
	# get the number of gene annotated with $f
	#
	my $s1 =  $self->getGOSingleCategoryNumbers($f);
	
	
	next if ($s1 == 0);
	next if ($self->{MAXCATEGORY} && ($s1 > $self->{MAXCATEGORY}));
	next if ($self->{MINCATEGORY} && ($s1 < $self->{MINCATEGORY}));

	
	my $p = Hypergeom::cumhyper($i,$s1,$s2,$self->{TOTALNBORFS});

	#print "Hypergeom::cumhyper($i,$s1,$s2,$self->{TOTALNBORFS}) $h_ref_names->{$f}->{NAME};\n";

	
	my %h_tmp = (
		     PVALUE  => $p, 
		     OVERLAP => $i, 
		     S1      => $s1, 
		     S2      => $s2, 
		     TEXT    => $h_ref_names->{$f}->{NAME}, 
		     NUM     => $f, 
		     N       => $self->{TOTALNBORFS},
		     EXP     => int(0.5 + $s2 * $s1 / $self->{TOTALNBORFS})
		     
		     );
	
	#
	# do nothing if p-vqlue is above threshold
	#
	#print "p=$p, $self->{PVALUE_THRESHOLD}, $self->{BONFERRONI}\n";
	next if (defined($self->{PVALUE_THRESHOLD}) && ( $p > $self->{PVALUE_THRESHOLD} )); 

	next if (($self->{BONFERRONI} == 1) && ($p * $self->{NBCATEGORIES} > 0.05));

	

	push @a_res, \%h_tmp;
    
	
    }
    
    my @a_res_bis = sort {$a->{PVALUE} <=> $b->{PVALUE}} @a_res;

    return \@a_res_bis;
}




1;
