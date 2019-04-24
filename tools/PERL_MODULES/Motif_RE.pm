#
#   handles motifs with regular expressions
#
use lib qw(/home/olly/PERL_MODULES);
package Motif_RE;
use Sets;

sub new {
    my ($self) = {};
    $self->{MOTIF}    = undef;
    $self->{MOTIFRE}  = undef;
    $self->{MOTIFTXT} = undef;
    $self->{SPECIES}  = [];
    $self->{FILE}     = undef;
    $self->{NAME}     = undef;

    bless($self);
    return $self;
}


sub loadTransfacMotif {
    my ($self, $f) = @_;

    $self->{FILE} = $f;

    $self->{NAME} = Sets::basename($f);

    
    open IN, $f or die "Pb openimg $f\n";
    my @a_motif = ();
    while (my $l = <IN>) {
	chomp $l;
	my $nn = substr($l, 0, 2);

	if ($nn =~ /BF/) {

	    my ($species) = $l =~ /Species\:\ (.+)$/;
	    push @{ $self->{SPECIES} }, $species;
	    
	}

	#  numbered line
	if ($nn =~ /\d\d/) {
	    my @a  = split /\s+/, $l;
	    push @a_motif, $a[5];
	}


	#print "$a[5]\n";
	
	
    }
    
    $self->{MOTIF}    =         \@a_motif;
    $self->{MOTIFTXT} = $self->getMotif;
    
    $self->{MOTIFRE} = $self->_motifToRE($self->{MOTIF});
    
    close IN;
}

#
#  create a motif using name and consensus
#
sub createMotifFromNameAndConsensus {
    my ($self, $n, $s) = @_;

    $self->{NAME}     = $n;
    $self->{MOTIF}    = [ split(//, $s) ];
    $self->{MOTIFTXT} = $s;        
    $self->{MOTIFRE}  = $self->_motifToRE($self->{MOTIF});
    
}


sub print {
    my ($self) = @_;
    
    print "$self->{NAME}\t";
    print $self->{MOTIFRE};
    print "\n";
}

sub hasSpecies {
    my ($self, $s) = @_;
    
    return scalar( grep (/$s/, @{ $self->{SPECIES} } ) );
}

sub getFileName {

    my ($self) = @_;

    if ($self->{FILE}) {
	return $self->{FILE};
    } else {
	return $self->{NAME};
    }
}


sub getName {

    my ($self) = @_;

    return $self->{NAME};
}


sub _motifToRE {
    my ($self, $a_ref) = @_;
    
    my $s_re = "";
    #my @a_letters = split //, $s_motif;

    foreach my $l (@$a_ref) {
        if    ($l =~ /[ATCG]/) { $s_re .= $l; }
        elsif ($l eq 'S')      { $s_re .= '[CG]'; }
        elsif ($l eq 'W')      { $s_re .= '[AT]'; }
        elsif ($l eq 'R')      { $s_re .= '[AG]'; }
        elsif ($l eq 'Y')      { $s_re .= '[CT]'; }
        elsif ($l eq 'K')      { $s_re .= '[GT]'; }
        elsif ($l eq 'M')      { $s_re .= '[AC]'; }
        elsif ($l eq 'B')      { $s_re .= '[CGT]'; }
        elsif ($l eq 'D')      { $s_re .= '[AGT]'; }
        elsif ($l eq 'H')      { $s_re .= '[ACT]'; }
        elsif ($l eq 'V')      { $s_re .= '[ACG]'; }
        elsif ($l eq 'N')      { $s_re .= '.'; }
    }
    
    

    return $s_re;

    
}


sub getRE {
    my ($self) = @_;
    
    return $self->{MOTIFRE}
    
    
}

sub getMotif {
    my ($self) = @_;
    
    return join("", @{$self->{MOTIF}})
}

#
# returns TRUE if given sequence matches the motif
#
sub match {
    
    my ($self, $s) = @_;
 
    # if the seq is smaller than the motif, embed it into a 
    #   much bigger one
    my $sc  = $s;
    my $scc = Sets::getComplement($sc);
    
    my $diff = length($self->{MOTIFTXT}) - length($sc);

#    print ">$self->{MOTIFTXT}\n";


    #
    #  the motif is larger than the sequence
    #
    
    #
    #  look at several substr of the motif
    #

    if ($diff > 0 ) {

	#$sc = ("N" x $diff) . $sc . ("N" x $diff);

	my $lenseq = length($sc);
	my $lenmot = length($self->{MOTIFTXT});

	my $l = $lenmot - $lenseq + 1;

	#print "l=$l\n";

	# slide the sequence along the motif
	for (my $i=0; $i<$l; $i++ ) {

	    # get the submotif starting at pos $i
	    my $subre = substr($self->{MOTIFTXT}, $i, $lenseq);
	    
	    # count the number of Ns in this substr
	    my $nbn = scalar(grep(/N/, split(//, $subre)));
	    
	    # this is the max nb of N for this seq
	    #my $max = 0.5*$lenseq+1;

	 #   print "Comparing $sc/$scc and $subre\n";


	    # if there are too many Ns, go on
	    #next if ($nbn > (0.5*$lenseq-1));
	    next if ($nbn > 1); #$lenseq-1);

	    #print "Comparing $sc/$scc and $subre\n";
	    
	    #print "Looking at subre $subre, nbn=$nbn ($max)\n";
	    my $re = $self->_motifToRE([ split(//, $subre) ]);
	    #print "Comparing $sc/$scc and $re\n";

	    #print "($re)\n";
	    if ( ($sc =~ /$re/) || ($scc =~ /$re/) ) {
		#print "MATCH\n";
		return 1;
	    } else {
		#print "\n";
		#return 0;
	    }
	    
	    
	}
    


    } else {



	my $re = $self->{MOTIFRE};
	#print "Comparing $sc/$scc and $re\n";

	if ( ($sc =~ /$re/) || ($scc =~ /$re/)) {
	    return 1;
	} else {
	    #return 0;
	}
    
    }
    
    return 0;
}



1;
