#
#   handles lists of motifs with regular expressions
#

package MotifLibrary;
use Sets;
use Motif;


sub new {
    my ($self) = {};
    $self->{OBJECTS} = [];
    $self->{VERBOSE} = 0;
    $self->{SPECIES} = [];
    $self->{MAXNBMOTIFS} = undef;
    bless($self);
    return $self;
}


sub setVerbose {
    my ($self, $l)   = @_;
    $self->{VERBOSE} = $l;
}

sub setMaxNbMotifs {
    my ($self, $l)   = @_;
    $self->{MAXNBMOTIFS} = $l;
}

#
#  load AlignACE motif output
#
sub loadAlignACEOutputFile {
    my ($self, $f) = @_;

    open ALI, $f or die "Could not open AlignACE file\n";
    
    while (1) {
	$l = <ALI>; last if (!$l);
	chomp $l;
	
	
	
	#print "$l\n";

	if ($l =~ /Input sequences\:/) {
	    
	    while (my $l ne "") {
		
		if (my ($ns, $ts) = $l =~ /\#(\d+)\t(.+?)$/) {

		    
		    $self->{SPECIES}->[$ns] = $ts;

		    #print "adding $ns $ts\n";
		}
		
		$l = <ALI>; chomp $l;
	    }
	    
	}


	if ($l =~ /Motif (\d+)/) {

	    $l = <ALI>; chomp $l;
	    
	    
	    my $mo = Motif->new;
	    
	    while (($l ne "") && ($l !~ /MAP Score\: [\d\.]+$/)) {

		if ($l =~ /\*/) {
		    $mo->setStars($l);
		} else {
		    my @a = split /\t/, $l;
		    
		    if (($a[0] !~ /[^N]N/) && ($a[0] !~ /N[^N]/)) {
		    
			$mo->addSiteFromSeq($a[0],  $a[1]);
			
		    }

		    #print "adding $a[0]\n";
		}

		$l = <ALI>; chomp $l;
	    }
	    
	    if (my ($sco) = $l =~ /MAP Score\: ([\d\.]+)$/) {
		$mo->setScore($sco);
	    }

	    

	    push @{  $self->{OBJECTS} }, $mo;

	    
	    last if (defined($self->{MAXNBMOTIFS}) && scalar(@{  $self->{OBJECTS} }) ==  $self->{MAXNBMOTIFS});

	}
	
    }

    close ALI;

}



sub sortByScore {
    my ($self) = @_;
    my @a = sort {$b->getScore() <=> $a->getScore()} @{  $self->{OBJECTS} };
    @{  $self->{OBJECTS} } = @a;
}






sub printInfo {
    my ($self) = @_;

    print "got " . scalar(@{  $self->{OBJECTS} }) . " motifs\n";

    print "score first motif = " . $self->{OBJECTS}->[0]->getScore() . "\n";
    
    
}

#
#  load a file full of AlignACE motifs
#
sub loadScanACEMotifs {
    my ($self, $l) = @_;
    
    my $a_ref = Sets::readSet($l);

    foreach my $s (@$a_ref) {
	my $m = Motif->new;
	$m->readScanACEMotif($s);
	$m->setName($s);
	push @{$self->{OBJECTS}}, $m;
    }

    if ($self->{VERBOSE} == 1) {
	print sprintf("Loaded %d motifs ..\n", scalar( @{$self->{OBJECTS}} ));
	print "\n";
    }
}


sub loadScanACEMotifList {
    my ($self, $l) = @_;
    
    my $a_ref = Sets::readSet($l);

    foreach my $s (@$a_ref) {
	my $m = Motif->new;
	$m->readScanACEMotif($s);
	$m->setName($s);
	push @{$self->{OBJECTS}}, $m;
    }

    if ($self->{VERBOSE} == 1) {
	print sprintf("Loaded %d motifs ..\n", scalar( @{$self->{OBJECTS}} ));
	print "\n";
    }
}





sub getMotifs {
    my ($self) = @_;

    return $self->{OBJECTS};
    
}


sub getOneMotif {
    my ($self, $n) = @_;
    
    return $self->{OBJECTS}->[$n];
}

#
#   match a given motif to the library
#   returns all matches
#
sub matchToAll {
    my ($self, $mym) = @_;
    
    if ($self->{VERBOSE} == 1) {
    
	print "Matches to :\n";
    }    

    my @a_matches = ();
    foreach my $m (@{$self->{OBJECTS}}) {
	my $score = $mym->compareACEScoreWithMotif($m);
	my @a_tmp = (Sets::basename($m->getName()), $score);

	if ($self->{VERBOSE} == 1) {
	    
	    print  $m->getName() . "\t" . $score . "\n"; 
	}
	push @a_matches, \@a_tmp;
    }

    
    my @a_sorted_matches = sort { $b->[1] <=> $a->[1] } @a_matches; 
    
    

    return \@a_sorted_matches;
}

1;
