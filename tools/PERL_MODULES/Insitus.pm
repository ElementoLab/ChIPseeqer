package Insitus;
use lib qw(/home/olly/PERL_MODULES);
use Database;
use Hypergeom;
use Sets;
use strict;


sub new {
    my $self  = {};
    

    $self->{DB} = Database->new;
    

    $self->{USER}             = "root";
    $self->{PASS}             = "";
    $self->{HOST}             = "";

    $self->{DB}->connect("FLYDB");

    $self->{ORF_SET}          = [];    

    $self->{TOTALNBORFS}      = 13522; #6200; #5538;  # Kellis et al number    

    $self->{PVALUE_THRESHOLD} = undef; 
    $self->{VERBOSE}          = 0;
    
    $self->{NBCATEGORIES}     = 0;
    
    $self->{BONFERRONI}       = 1;



    bless($self); 
    return $self;
    
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



#
# get all annotations for a given gene
#
sub getGeneAnnotations {

    my ($self, $orf) = @_;
	
    my $sql   = undef;
    
    $sql = "SELECT * FROM INSITU_ANNOTATIONS WHERE ID  = '$orf'";
	
    my $a_ref = $self->{DB}->queryAllRecordsRef($sql);
    
    return $a_ref;
    
}





#
#  get all info related to a given tissue id : text and count
#
sub getInfoTissues {
    
    my ($self) = @_;
    
    my $sql1          = "select TISSUEID, count(*) as COUNT from INSITU_ANNOTATIONS group by TISSUEID";
    my $a_ref1        = $self->{DB}->queryAllRecordsRef($sql1);
    my $h_ref_tissues = Sets::SQLRefToIndex($a_ref1, "TISSUEID");

    my $sql2          = "select * from INSITU_TISSUES";
    my $a_ref2        = $self->{DB}->queryAllRecordsRef($sql2);
    foreach my $r (@$a_ref2) {
	$h_ref_tissues->{$r->{TISSUEID}}->{TISSUE} = $r->{TISSUE};
    }
    
    return $h_ref_tissues;
    
}


#
#  get all info related to a given stage id : text and count
#
sub getInfoStages {
    
    my ($self) = @_;
    
    my $sql1          = "select distinct ID, STAGEID from INSITU_ANNOTATIONS;";
    
    #print "$sql1\n";
    my $a_ref1        = $self->{DB}->queryAllRecordsRef($sql1);
    my %h = ();
    my $h_ref_stages = \%h;
    foreach my $r (@$a_ref1) {
	$h_ref_stages->{ $r->{STAGEID} }->{COUNT}  ++;   #Sets::SQLRefToIndex($a_ref1, "STAGEID");
    }

    $self->{NBSTAGES} = scalar(keys(%$h_ref_stages));

    #foreach my $k (keys(%$h_ref_stages)) {
    #	print "$k\t$h_ref_stages->{$k}->{COUNT}\n";
    #}
    
    #print "toto\n";

    my $sql2          = "select * from INSITU_STAGES";
    my $a_ref2        = $self->{DB}->queryAllRecordsRef($sql2);
    foreach my $r (@$a_ref2) {
	$h_ref_stages->{$r->{STAGEID}}->{STAGE} = $r->{STAGE};
    }
    
    return $h_ref_stages;
    
}



#
#  get all info related to tissues : text and count
#
sub getInfoTissues {
    
    my ($self) = @_;
    
    my $sql1          = "select distinct ID, TISSUEID from INSITU_ANNOTATIONS;";
    
    #print "$sql1\n";
    my $a_ref1        = $self->{DB}->queryAllRecordsRef($sql1);
    my %h = ();
    my $h_ref_tissues = \%h;
    foreach my $r (@$a_ref1) {
	$h_ref_tissues->{ $r->{TISSUEID} }->{COUNT}  ++;   #Sets::SQLRefToIndex($a_ref1, "STAGEID");
    }

    $self->{NBTISSUES} = scalar(keys(%$h_ref_tissues));
    
    #foreach my $k (keys(%$h_ref_tissues)) {
    #print "$k\t$h_ref_tissues->{$k}->{COUNT}\n";
    #}
    
    #print "Tissues\n";

    my $sql2          = "select * from INSITU_TISSUES";
    my $a_ref2        = $self->{DB}->queryAllRecordsRef($sql2);
    foreach my $r (@$a_ref2) {
	$h_ref_tissues->{$r->{TISSUEID}}->{TISSUE} = $r->{TISSUE};
    }
    
    return $h_ref_tissues;
    
}


sub getStagesAndTissues {
    
    my ($self) = @_;

    #  get counts
    my $sql1          = "select distinct INSITU_ANNOTATIONS.STAGEID, STAGE, INSITU_ANNOTATIONS.TISSUEID, TISSUE from INSITU_ANNOTATIONS, INSITU_TISSUES, INSITU_STAGES where INSITU_TISSUES.TISSUEID = INSITU_ANNOTATIONS.TISSUEID and INSITU_STAGES.STAGEID = INSITU_ANNOTATIONS.STAGEID";
    my $a_ref1        = $self->{DB}->queryAllRecordsRef($sql1);
    
    return $a_ref1;

 
    
}


sub getGenesAnnotatedWithStageTissue {

     my ($self, $s, $t) = @_;
     
     
     my $sql1          = "select distinct ID from INSITU_ANNOTATIONS where TISSUEID = '$t' and STAGEID = '$s'"; 
     my $a_ref1        = $self->{DB}->queryAllRecordsRef($sql1);
     return $a_ref1;
     

}


sub getInfoStagesAndTissues {
    
    my ($self) = @_;

    #  get counts
    my $sql1          = "select distinct ID, STAGEID, TISSUEID from INSITU_ANNOTATIONS";
    my $a_ref1        = $self->{DB}->queryAllRecordsRef($sql1);
    
    my %h = ();
    my $h_ref_stages_tissues = \%h;

    #  get tissue names 
    my $sql2          = "select * from INSITU_TISSUES";
    my $a_ref2        = $self->{DB}->queryAllRecordsRef($sql2);
    my $h_ref2        = Sets::SQLRefToIndex($a_ref2, "TISSUEID");
    
    #  get stage names
    my $sql3          = "select * from INSITU_STAGES";
    my $a_ref3        = $self->{DB}->queryAllRecordsRef($sql3);
    my $h_ref3        = Sets::SQLRefToIndex($a_ref3, "STAGEID");

    foreach my $r (@$a_ref1) {
	my $k = $r->{STAGEID} . "_" . $r->{TISSUEID};
	$h_ref_stages_tissues->{ $k }->{COUNT} ++;
	$h_ref_stages_tissues->{ $k }->{TISSUE} = $h_ref2->{ $r->{TISSUEID} }->{TISSUE};
	$h_ref_stages_tissues->{ $k }->{STAGE } = $h_ref3->{ $r->{STAGEID } }->{STAGE };
    }
    
    $self->{NBSTAGESTISSUES} = scalar(keys(%$h_ref_stages_tissues));
    

    return $h_ref_stages_tissues;
    
}



#   
#  get the STAGE enrichment of an array of ORFS
#
sub calculateEnrichments {
    
    my ($self) = @_;
    
    if (scalar(@{$self->{ORF_SET}}) == 0) {
	die "Please provide a non-empty array\n";
    }
    

    #  
    #   get info for all tissues, and all stages and all combinations of stages/tissues
    #
    my $h_ref_tissues = $self->getInfoTissues();
    my $h_ref_stages  = $self->getInfoStages();
    my $h_ref_stages_tissues  = $self->getInfoStagesAndTissues();
    
    #
    #   for each ORF in the set, get the function with which it is annotated 
    #
    my %stage_count  = ();
    my %tissue_count  = ();
    my %stage_tissue_count  = ();

    foreach my $o (@{$self->{ORF_SET}}) {

	#print "Gene $o\n" if ($self->{VERBOSE});
	#
	# get all the categories
	#
	my $a_ref_annot  = $self->getGeneAnnotations($o);

	#
	# get the different stages in which this gene is expressed
	#
	
	my %stage_there = ();
	my %tissue_there = ();
	my %stage_tissue_there = ();
	
		
	foreach my $r (@$a_ref_annot) {
	    
	    $stage_there{ $r->{STAGEID} } = 1;
	    $tissue_there{ $r->{TISSUEID} } = 1;
	    $stage_tissue_there{ $r->{STAGEID} . "_" . $r->{TISSUEID} } = 1;
	    
	    #$tissue_count{ $r->{TISSUEID} } ++;
	    #$stage_count { $r->{STAGEID } } ++;
	    #$stage_tissue_count { $r->{STAGEID } }{ $r->{TISSUEID } } ++;
	    
	    
	    
	}
	
	foreach my $stageid (keys(%stage_there)) {
	    $stage_count { $stageid } ++;
	}

	foreach my $tissueid (keys(%tissue_there)) {
	    $tissue_count { $tissueid } ++;
	}

	foreach my $stage_tissue_id (keys(%stage_tissue_there)) {
	    $stage_tissue_count { $stage_tissue_id } ++;
	}

	#$i_cnt++;
    }

    #
    # now traverse the non-zero categories and calc a p-value
    #
    #my @a_res = ();


    my $s2 = scalar(@{$self->{ORF_SET}});

    #
    #  stages
    #
    $self->{STAGES} = [];
    while (my ($stageid,$count) = each(%stage_count)) {

	#print "$stageid,$count\n";
    	
	#
	# get the total number of genes annotated with this $stageid
	#
	my $s1 =  $h_ref_stages->{$stageid}->{COUNT};
	
	
	next if ($s1 == 0);
		
	my $p = Hypergeom::cumhyper($count,$s1,$s2,$self->{TOTALNBORFS});
	
	my %h_tmp = (
		     PVALUE  => $p, 
		     OVERLAP => $count, 
		     S1      => $s1, 
		     S2      => $s2, 
		     STAGE   => $h_ref_stages->{$stageid}->{STAGE}, 
		     NUM     => $stageid, 
		     N       => $self->{TOTALNBORFS},
		     EXP     => int(0.5 + $s2 * $s1 / $self->{TOTALNBORFS})
		     );
	
#	print "$p\n";

	#
	# do nothing if p-vqlue is above threshold
	#
	next if (defined($self->{PVALUE_THRESHOLD}) && ( $p > $self->{PVALUE_THRESHOLD} )); 

	next if (defined($self->{BONFERRONI}) && ($p * $self->{NBSTAGES} > 0.05));

	push @{ $self->{STAGES} }, \%h_tmp;
    
    }
    
    @{ $self->{STAGES} } = sort {$a->{PVALUE} <=> $b->{PVALUE}} @{ $self->{STAGES} };





    #
    #  TISSUES
    #
    $self->{TISSUES} = [];
    while (my ($tissueid,$count) = each(%tissue_count)) {
    	
	#
	# get the total number of genes annotated with this $stageid
	#
	my $s1 =  $h_ref_tissues->{$tissueid}->{COUNT};
	
	
	next if ($s1 == 0);
		
	my $p = Hypergeom::cumhyper($count,$s1,$s2,$self->{TOTALNBORFS});
	
	my %h_tmp = (
		     PVALUE  => $p, 
		     OVERLAP => $count, 
		     S1      => $s1, 
		     S2      => $s2, 
		     TISSUE  => $h_ref_tissues->{$tissueid}->{TISSUE}, 
		     NUM     => $tissueid, 
		     N       => $self->{TOTALNBORFS}
		     );
	
	#
	# do nothing if p-vqlue is above threshold
	#
	next if (defined($self->{PVALUE_THRESHOLD}) && ( $p > $self->{PVALUE_THRESHOLD} )); 

	next if (defined($self->{BONFERRONI}) && ($p * $self->{NBSTAGESTISSUES} > 0.05));

	push @{ $self->{TISSUES} }, \%h_tmp;
    
    }
    
    @{ $self->{TISSUES} } = sort {$a->{PVALUE} <=> $b->{PVALUE}} @{ $self->{TISSUES} };


    
    
    #
    #  STAGE / TISSUES
    #
    $self->{STAGES_TISSUES} = [];
    while (my ($stage_tissue_id,$count) = each(%stage_tissue_count)) {
    	
	#
	# get the total number of genes annotated with this $stageid
	#
	my $s1 =  $h_ref_stages_tissues->{$stage_tissue_id}->{COUNT};
	
	
	next if ($s1 == 0);
		
	my $p = Hypergeom::cumhyper($count,$s1,$s2,$self->{TOTALNBORFS});
	
	my %h_tmp = (
		     PVALUE  => $p, 
		     OVERLAP => $count, 
		     S1      => $s1, 
		     S2      => $s2, 
		     STAGE_TISSUE  => $h_ref_stages_tissues->{$stage_tissue_id}->{STAGE} . ":" . $h_ref_stages_tissues->{$stage_tissue_id}->{TISSUE}, 
		     NUM     => $stage_tissue_id, 
		     N       => $self->{TOTALNBORFS},
		     EXP     => int(0.5 + $s2 * $s1 / $self->{TOTALNBORFS})

		     );
	
	#
	# do nothing if p-vqlue is above threshold
	#
	next if (defined($self->{PVALUE_THRESHOLD}) && ( $p > $self->{PVALUE_THRESHOLD} )); 

	next if (defined($self->{BONFERRONI}) && ($p * $self->{NBTISSUES} > 0.05));

	push @{ $self->{STAGES_TISSUES} }, \%h_tmp;
    
    }
    
    @{ $self->{STAGES_TISSUES} } = sort {$a->{PVALUE} <=> $b->{PVALUE}} @{ $self->{STAGES_TISSUES} };




}

sub getStageTissueEnrichments {
    my ($self) = @_;
    
    return $self->{STAGES_TISSUES};
}

sub getStageEnrichments {
    my ($self) = @_;
    
    return $self->{STAGES};
}


sub getTissueEnrichments {
    my ($self) = @_;
    
    return $self->{TISSUES};
}


1;
