#
#   handles lists of motifs with regular expressions
#

package Library_Motif_RE;
use Sets;
use Motif_RE;
use Table;


sub new {
    my ($self) = {};
    $self->{OBJECTS} = [];
    $self->{SPECIES} = undef;
    $self->{VERBOSE} = 0;
    my %IDX = ();
    $self->{OCCURRENCES} = \%IDX;

    bless($self);
    return $self;
}

sub setSpecies {
    
    my ($self, $l) = @_;
    $self->{SPECIES} = $l;
    
}

sub loadMotifList {
    my ($self, $l) = @_;
    
    my $a_ref = Sets::readSet($l);

    foreach my $s (@$a_ref) {
	my $m = Motif_RE->new;
	$m->loadTransfacMotif($s);
	
	if ($self->{SPECIES} && $m->hasSpecies($self->{SPECIES})) {
	    push @{$self->{OBJECTS}}, $m;
	} else {
	    push @{$self->{OBJECTS}}, $m;
	}
    }

    if ($self->{VERBOSE} == 1) {
	print sprintf("Loaded %d motifs ..\n", scalar( @{$self->{OBJECTS}} ));
	print "\n";
    }
}

#
# en gros, charge la liste Kellis et al.
#
sub loadSimpleMotifList {
    my ($self, $f) = @_;

    my $ta = Table->new;
    $ta->loadFile($f);

    my $a_ref = $ta->getArray();
    foreach my $r (@$a_ref) {
	
	# create a new motif
	my $m = Motif_RE->new;
	
	$m->createMotifFromNameAndConsensus($r->[0], $r->[1]);
	
	#$m->print();
    
	push @{$self->{OBJECTS}}, $m;
    }

}


#
#  load a table name => a_ref (name, nb, z)
#
sub loadMotifOccurences {
    my ($self, $l) = @_;
    
    my $ta = Table->new;
    $ta->loadFile($l);

    $self->{OCCURRENCES} = $ta->getIndex(0);
    
    
}


sub getNbOccurrences {
    my ($self, $m) = @_;

    return $self->{OCCURRENCES}->{$m}->[1];
    
}


sub getMotifs {
    my ($self) = @_;

    return $self->{OBJECTS};
    
}

sub matchToAll {
    my ($self, $s) = @_;
    
    my @a_matches = ();
    foreach my $m (@{$self->{OBJECTS}}) {
	   #x print $m->getFileName() . "\n";

	if ($m->match($s)) {
	    #print "V: " . $m->getName() . "->" . $self->getNbOccurrences($m->getName()) . "\n";
	    my @a_tmp = ($m->getFileName(), $m->getMotif(), $self->getNbOccurrences($m->getName()));
	    push @a_matches,  \@a_tmp;
	}
	
    }
    
    return \@a_matches;
}

1;
