package Yeast;

use lib qw(/home/olly/PERL_MODULES);
use Hypergeom;

use Database;

sub new {    
    my $self  = {};
 
    $self->{DB} = Database->new;
    $self->{DB}->connect("YEASTHOUSE") or die "Cannot connect\n";
    #$self->{DB}->setVerbose(1);
    $self->{VERBOSE} = 0;
    $self->{TOTAL_NB_ORFS} = undef;
    bless($self);           # but see below
    return $self;    

}

sub setTotalNbORFs {
    my ($self, $o) = @_;
    $self->{TOTAL_NB_ORFS} = $o;
}

sub setVerbose {
    my ($self, $o) = @_;
    $self->{VERBOSE} = $o;

    $self->{DB}->setVerbose($o);
    

}


sub getGeneName {

    my ($self, $o) = @_;
    

    $s = "select * from YEASTHOUSE where ORF='$o'";
    
    
    my $a = $self->{DB}->queryOneRecord($s);
    
    if ($a->{GENENAME} =~ /ORF/) {
	return $o;
    } else {
	return  $a->{GENENAME};
    }
}


sub getORFName {

    my ($self, $g) = @_;
    

    $s = "select * from YEASTHOUSE where GENENAME='$g'";
    
    
    my $a = $self->{DB}->queryOneRecord($s);
    
    return  ($a->{ORF}?$a->{ORF}:$g);

}



sub getORF {

    my ($self, $o) = @_;
    

    $s = "select * from YEASTHOUSE where ORF='$o'";
    
    
    my $a = $self->{DB}->queryOneRecord($s);
    
    return  $a;

}


sub getAllProcesses {
    
    my ($self, $o) = @_;
      
    my $s = "select GO_ANNOTATION.* from GO_ANNOTATION where GO_ANNOTATION.ORF='$o' and TYPE='P'";
    
    #print $s;
    
    my @a = $self->{DB}->queryAllRecords($s);
        
    return  \@a;
}





sub getGeneticInteractions {

    my ($self, $o) = @_;
    
    
    my $s = "select  *  from PROTEINDNA_INTERACTIONS where P1 = '$o'";
    
    #print $s;
    
    my @a = $self->{DB}->queryAllRecords($s);
    
    return \@a;
}


sub getAllInteractionsByType {

    my ($self, $exp, $src) = @_;
    
    
    my $s = "select distinctrow P1,P2  from INTERACTIONS where EXP = '$exp'";
    
    if ($src) {
	$s .= " and SRC = '$src'";
    }
    
    #print $s;
    
    my @a = $self->{DB}->queryAllRecords($s);
    
    return \@a;
}


sub getPhysicalInteractions {

    my ($self, $o, $src) = @_;
    
    
    my $s = "select distinct P2 from INTERACTIONS where P1 = '$o' AND EXP != 'Synthetic Lethality'";
    
    if ($src) {
	$s .= " and SRC = '$src'";
    }
    
    #print $s;
    
    my $a = $self->{DB}->queryAllRecordsRef($s);
    
    my $b = Sets::SQLRefToSet($a, "P2");

    return $b;
}


sub getPhenotypes {
    my ($self, $o) = @_;

    my $ca = Categories->new;
    $ca->setResource("/home/olly/DATA/YEASTS/PHENOTYPE");
    
    
    
    return $ca->getCategoriesText($o);

}


#
# get all TF
#
sub getAllTFs {

    my ($self) = @_;
    
    #my $s = "select distinct YEASTHOUSE.GENENAME, YEASTHOUSE.ORF  from GO_ANNOTATION, YEASTHOUSE where ((GOID = 30528) OR (GOID = 3677)) AND YEASTHOUSE.ORF=GO_ANNOTATION.ORF order by YEASTHOUSE.GENENAME";
    my $s = "select  *  from FACTORS";
    
   # print $s;
    
    my @a = $self->{DB}->queryAllRecords($s);
    
    


    return  \@a;
}

#
# get the correlation value between two gene names/orfs
#


#
#  get divergence for all genes 
#

sub getDivergences {
    my ($self) = @_;
    
    my $s1 = "select distinct ORF, DISTANCE_JC_BAYANUS from YEASTHOUSE where DISTANCE_JC_BAYANUS > 0.0";
    
    
    my $a_ref = $self->{DB}->queryAllRecordsRef($s1);
   
    return $a_ref;
    
}


#
#  get divergence for all genes 
#

sub getYeastFeature {
    my ($self, $f) = @_;
    
    my $s1 = "select distinct ORF, $f from YEASTHOUSE where $f > 0.0";
    
    
    my $a_ref = $self->{DB}->queryAllRecordsRef($s1);
   
    return $a_ref;
    
}


#
#  get essential genes 
#
sub getEssentialGenes {
    my ($self) = @_;
    
    my $s1 = "select distinct ORF from YEASTHOUSE where ESSENTIAL = 1";
    
    
    my $a_ess = $self->{DB}->queryAllRecordsRef($s1);
   
    my $a_set = Sets::SQLRefToSet($a_ess, "ORF");

    return $a_set;
    
}


sub getRegulators {
    my ($self, $o, $c, $s, $p) = @_;

    my $sql = "select distinctrow REGULATOR  from CHIP203 where ORF = '$o'";

    if ($c) {
	$sql .= " and CONDITION = '$c';";
    }

    #if ($c) {
	$sql .= " and PVALUE <= 0.001;";
    #}
    
    my $a_reg = $self->{DB}->queryAllRecordsRef($sql);
    
    my $a_ref = Sets::SQLRefToSet($a_reg, "REGULATOR");

    return Sets::removeDuplicates($a_ref);
    
}


sub getNumberOfTargets {
    my ($self, $t) = @_;

    my $sql = "select count(*) as COUNT  from CHIP203 where REGULATOR = '$t'";
    
    my $a = $self->{DB}->queryOneRecord($sql);
    
    return $a->{COUNT};
}

#
#  get the CHIP enrichment of a group of genes
#  
#
sub getCHIPEnrichment {
    my ($self, $a_ref_set, $totalnbgenes, $pv) = @_;
    
    # get all ChIP regulators
    #my $sql = "select distinct REGULATOR  from CHIP203";
    my $sql = "select distinct REGULATOR  from CHIP";
        
    my $a_reg = $self->{DB}->queryAllRecordsRef($sql);
    my @a_res = ();

    foreach my $r (@$a_reg) {
	
	# for each regulator, get the set of regulated genes
	#$sql          = "select ORF from CHIP203 where REGULATOR = '$r->{REGULATOR}' and PVALUE < 0.001"; 
	$sql          = "select ORF from CHIP where REGULATOR = '$r->{REGULATOR}' and PVALUE < 0.001"; 

	my $a_orfs_db = $self->{DB}->queryAllRecordsRef($sql);

	my @a_orfs = ();
	foreach my $d_ref (@$a_orfs_db) {
	    push @a_orfs, $d_ref->{ORF};
	}

	# calculate the overlap and its significance
	my $a_ovl = Sets::getOverlapSet($a_ref_set, \@a_orfs);

	my $s1    = scalar(@$a_ref_set);
	my $s2    = scalar(@a_orfs);
	my $i     = scalar(@$a_ovl);
	
	my $p     = Hypergeom::cumhyper($i, $s1, $s2, $totalnbgenes);

	next if ($p > $pv);


	my %h_tmp = (REGULATOR => $r->{REGULATOR}, 
		     P         => $p,
		     S1        => $s1,
		     S2        => $s2);

	push @a_res, \%h_tmp;
    }

    
    my @a_res_sorted = sort { $a->{P} <=> $b->{P} } @a_res;


    # returns the results
    return \@a_res_sorted;
}

#
#  kmer = given kmer to be transformed
#  cpt  = compareace threshold 
#
sub getMotifEnrichmentAndCompareACEscoreForKmer {
    my ($self, $a_ref, $totalnbgenes, $pv, $kmer, $cpt) = @_;
	
    my $a_res = $self->getMotifEnrichment($a_ref, $totalnbgenes, $pv);

    # for each motif, get the compareace score with the kmer

    my $mo2 = Motif->new;
    $mo2->createMotifFromSite($kmer, '*' * length($kmer));

    my @out = ();
    foreach my $r (@$a_res) {
	my $f = "/home/olly/DATA/YEASTS/KNOWN_MOTIFS/$r->[0].ace";
	
	#print "getting $f\n";

	my $mo1 = Motif->new;
	$mo1->readScanACEMotif($f);

	my $sc = $mo1->compareACEScoreWithMotif($mo2);

	#print "sc=$sc\n";
	
	if ($sc >= $cpt) {
	    my @a = (@$r, $sc);

	    push @out, \@a;
	}
	
    }


    return \@out;
}

#
#  get motif enrichment
#
sub getMotifEnrichment {
    my ($self, $a_ref, $totalnbgenes, $pv) = @_;
    
    # get the pre-computed genes
    my $a_files = Sets::getFiles("/home/olly/PROGRAMS/MOTIF_TARGETS/FILES/*.txt");

    

    my @a_res = ();
    
    foreach my $f (@$a_files) {
	#print "$f\n";
	my $a_ref_set = Sets::readSet($f);
	my $a_ovl     = Sets::getOverlapSet($a_ref_set, $a_ref);
	my $s1               = scalar(@$a_ref_set);
	my $s2               = scalar(@$a_ref);
	my $i                = scalar(@$a_ovl);
	my $p                = Hypergeom::cumhyper($i, $s1, $s2, $totalnbgenes);
	
	if ($p < $pv) {
	    my @a_tmp = (Sets::basename($f), $p);
	    push @a_res, \@a_tmp;
	    #print sprintf("%s\t%e\n", Sets::basename($f), $p);
	} 
    }

    my @a_res_sorted = sort { $a->[1] <=> $b->[1] } @a_res;

    return \@a_res_sorted;

}










sub getCHIPdata {
    my ($self, $r, $o) = @_;

    my $s = "select  *  from CHIP where REGULATOR='$r' and ORF='$o'";
    
    #print $s;
    
    my $a = $self->{DB}->queryOneRecord($s);
    
    


    return  $a;
}

#
# get all the ORFs regulated by a given TF above a given ratio
#
sub getAllCHIP_above_ratio {
    
    my ($self, $tf, $ratio) = @_;

    my $s = "select  *  from CHIP where REGULATOR='$tf' and RATIO>='$ratio'";
 
    #print "$s\n";
   
    my @a = $self->{DB}->queryAllRecords($s);

    return  \@a;
    
}




sub getCorrelationData {
    my ($self, $o) = @_;

    my $s = "select YEASTHOUSE.GENENAME, YEASTHOUSE.ORF, CORRELATION   from CORRELATIONS, YEASTHOUSE where REGULATOR='$o' and YEASTHOUSE.ORF = CORRELATIONS.ORF order by CORRELATION desc";
    
    #print $s;
    
    my @a = $self->{DB}->queryAllRecords($s);
    
    return  \@a;
}


sub getMotifs {
    my ($self, $r, $o) = @_;

    my $s = "select *   from MOTIFS where REGULATOR='$r' and ORF = '$o'";
    
    #print $s;
    
    my @a = $self->{DB}->queryAllRecords($s);
    
    return  \@a;
}






sub getCorrelation {

    #my ($self, $g1, ) = @_;

}



1;
