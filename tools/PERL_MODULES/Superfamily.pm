package Superfamily;
use Database;


sub new {
    my $self  = {};
    $self->{DB}      = Database->new;
    $self->{USER}    = "root";
    $self->{PASS}    = "";
    $self->{HOST}    = "localhost";
    
    $self->{CONNECTED} = 0;

    $self->{DB}->setID($self->{USER}, $self->{PASS}, $self->{HOST});
    
    bless($self); 
    return $self;
    
    
}

sub _connect {
    my ($self) = @_;

    if ($self->{CONNECTED} == 0) {
	$self->{DB}->connect("SUPERFAMILY");
    }
    
}

sub setID {
    my ($self, $u, $p, $h) = @_;
    
    $self->{USER} = $u;
    $self->{PASS} = $p;
    $self->{HOST} = $h;

    $self->{DB}->setID($self->{USER}, $self->{PASS}, $self->{HOST});

    
}

#
#  returns an array of domains for a given protein
#
sub getDomains {
    
    my ($self, $o) = @_;
    
    $self->_connect;
    
    my $s = "select ASSIGNMENTS.*, MODELS.TEXT from ASSIGNMENTS, MODELS where SEQUENCEID = '$o' and MODELS.MODELID = ASSIGNMENTS.MODELID";
    
    my $a_ref = $self->{DB}->queryAllRecordsRef($s);
    
    return $a_ref;
}


sub  getSFAMIDs {
    my ($self, $a_ref) = @_;
    
    my @a_set = ();
    foreach my $r (@$a_ref) {
	push @a_set, $r->{SFAMID};
    }

    return \@a_set;

}


#
#  delete tandem domains
#
sub  getStackOfSFAMIDs {
    my ($self, $a_ref) = @_;
    
    my $a_ref_dom = $self->getSFAMIDs($a_ref);
    
    my @a_stack = ();
    foreach my $r (@$a_ref_dom) {
	push @a_stack, $r if ($a_stack[$#a_stack] != $r);
    }

    return \@a_stack;


}


sub  sortDomainsByRegion {

    my ($self, $a_ref) = @_;
    
    my @a_new = sort _sortByRegion @$a_ref;

    return \@a_new;
}

sub  _sortByRegion {

    my $ar  = $a->{REGION};
    my $br  = $b->{REGION};

    my ($ars) = $ar =~ /^(\d+)\-/;
    my ($brs) = $br =~ /^(\d+)\-/;

    #print "$ars\t$brs\n";

    return $ars <=> $brs; 
}



#
#  get all the genomes
#
sub getAllGENOMEID {

    my ($self) = @_;
    
    $self->_connect;
    
    my $s = "select distinct GENOMEID from ASSIGNMENTS";
    
    my $a_ref = $self->{DB}->queryAllRecordsRef($s);
    
    return $a_ref;
    
}


#
#  get all the SFAMID
#
sub getAllSFAMID {

    my ($self) = @_;

    
    $self->_connect;
    
    
    my $s = "select distinct SFAMID from ASSIGNMENTS";
    
    my $a_ref = $self->{DB}->queryAllRecordsRef($s);
    
    return $a_ref;
    
}

sub getNbAssignments {
    
    my ($self, $s) = @_;
	

    $self->_connect;
    
    
    my $s = "select GENOMEID, count(*) as COUNT from ASSIGNMENTS group by GENOMEID";
    
    my $a_ref = $self->{DB}->queryAllRecordsRef($s);
    
    return $a_ref;
}

#
#   get the number of genes associated with 1 SFAMID / genome
#



#
#   get the number of occurences of a SFAMID / genome
#
sub getNbSFAMIDperGenome {

    my ($self, $d) = @_;
	

    $self->_connect;
    
    
    my $s = "select GENOMEID,count(*) as COUNT from ASSIGNMENTS where SFAMID = '$d' group by GENOMEID";
    
    my $a_ref = $self->{DB}->queryAllRecordsRef($s);
    
    return $a_ref;
    
}


sub getAllGenesBySFAMID {
    my ($self, $d, $o) = @_;

    $self->_connect;
    
    
    my $s = "select SEQUENCEID from ASSIGNMENTS where SFAMID = '$d'";
   
    if (defined($o)) {
	$s .= " and GENOMEID = '$o'";
    }
 
    print $s;

    my $a_ref = $self->{DB}->queryAllRecordsRef($s);
    
    return $a_ref;
    

}

1;
