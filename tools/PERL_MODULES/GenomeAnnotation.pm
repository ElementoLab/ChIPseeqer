package GenomeAnnotation;
use Table;
use Sets;
use strict;

sub new {
    my ($self)  = {};
    $self->{FILE} = undef;
 
       
    $self->{GENENAMES}      = [];
    $self->{CHR}            = Sets::initHash;
    $self->{EXONS}          = Sets::initHash;
    $self->{BOUNDARIES}     = Sets::initHash;
    $self->{BOUNDARIES_CDS} = Sets::initHash;
    $self->{STRAND}         = Sets::initHash;


    bless($self);
    return $self;

}


sub read {
    my ($self, $s) = @_;
    
    my $ta = Table->new;

    $ta->loadFile($s);
    
    my $a_ref = $ta->getArray();
 
    
    foreach my $r (@$a_ref) {
	push @{ $self->{GENENAMES} }, $r->[0];

	my @a_exon = ($r->[2], $r->[3]);
	push @{ $self->{EXONS}->{ $r->[0] } }, \@a_exon;
	$self->{CHR}->{ $r->[0] }             = $r->[1]; 
	$self->{BOUNDARIES}->{ $r->[0] }->[0] = (defined($self->{BOUNDARIES}->{ $r->[0] }->[0])?Sets::min($self->{BOUNDARIES}->{ $r->[0] }->[0], $r->[2]):$r->[2]);
	$self->{BOUNDARIES}->{ $r->[0] }->[1] = (defined($self->{BOUNDARIES}->{ $r->[0] }->[1])?Sets::max($self->{BOUNDARIES}->{ $r->[0] }->[1], $r->[3]):$r->[3]);
	$self->{STRAND}->{ $r->[0] }          = $r->[4];
	
	# no need to do that if this exon is not part of the CDS
	next if (($r->[5] eq "") && ($r->[6] eq ""));
		
	$self->{BOUNDARIES_CDS}->{ $r->[0] }->[0] = (defined($self->{BOUNDARIES_CDS}->{ $r->[0] }->[0])?Sets::min($self->{BOUNDARIES_CDS}->{ $r->[0] }->[0], $r->[5]):$r->[5]);
	$self->{BOUNDARIES_CDS}->{ $r->[0] }->[1] = (defined($self->{BOUNDARIES_CDS}->{ $r->[0] }->[1])?Sets::max($self->{BOUNDARIES_CDS}->{ $r->[0] }->[1], $r->[6]):$r->[6]);

	#print "$r->[0] CDS= " . $self->{BOUNDARIES_CDS}->{ $r->[0] }->[0] . " / " . $self->{BOUNDARIES_CDS}->{ $r->[0] }->[1] . "\n";
    }


    #
    #  correct the CDS boundaries if not annotated
    #
    
    foreach my $r (@{ $self->{GENENAMES} } ) {
	if (!defined($self->{BOUNDARIES_CDS}->{ $r })) {
	    $self->{BOUNDARIES_CDS}->{ $r }->[0] = $self->{BOUNDARIES}->{ $r }->[0];
	    $self->{BOUNDARIES_CDS}->{ $r }->[1] = $self->{BOUNDARIES}->{ $r }->[1];
	}
    }
}

sub getStrand {
    my ($self, $g) = @_;

    return $self->{STRAND}->{ $g };
    
    
}


sub getBoundariesCDS {
    my ($self, $g) = @_;

    return $self->{BOUNDARIES_CDS}->{ $g };
    
    
}

#
# takes a list of genes, sort them by chr, order them by their start
#
#
sub orderGenesOnChromosomes {
    my ($self, $a_ref) = @_;

    my %H = ();
    foreach my $r (@$a_ref) {
	push @{ $H{ $self->{CHR}->{ $r } } }, $r; 
    } 
    
    #print "ordering genes\n";

    foreach my $k (keys(%H)) {
	
	#print scalar(@{ $H{ $k } }); print " genes on $k \n";

	my @t = ();
	foreach my $g (@{ $H{ $k } }) {
	    my @a_tmp = ($g, $self->{BOUNDARIES_CDS}->{$g}->[0]);
	    
	    #print "x for '$g'=" . $self->{BOUNDARIES_CDS}->{$g}->[0] . "\n";

	    push @t, \@a_tmp;
	}

	@t = sort { $a->[1] <=> $b->[1] } @t;


	@{ $H{ $k } } = ();
	foreach my $g (@t) {
	    push @{ $H{ $k } }, $g->[0];
	}

	

    }

    return \%H;
    
}


sub genomicDistanceBetweenGenes {
    
    my ($self, $g1, $g2) = @_;
    
    my $c1 = $self->{CHR}->{$g1};
    my $c2 = $self->{CHR}->{$g2};

    if ($c1 != $c2) {
	return undef;
    }

    my $b1 = $self->{BOUNDARIES_CDS}->{$g1};
    my $b2 = $self->{BOUNDARIES_CDS}->{$g2};

    #print "$b1->[0], $b1->[1], $b2->[0], $b2->[1]\n";

    if (Sets::sequencesOverlap($b1->[0], $b1->[1], $b2->[0], $b2->[1])) {
	return -1;
    }
    
    # which one if first 
    if ($b1->[1] < $b2->[0]) {
	return $b2->[0] - $b1->[1];
    } else {
	return $b1->[0] - $b2->[1];
    }

}


1;
