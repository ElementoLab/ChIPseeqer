package Vienna;
use Sets;

sub new {
    my ($class) = @_;
    my ($self) = {};

    $self->{RNAFOLD} = "/home/elemento/PERL_MODULES/PROGRAMS/ViennaRNA-1.6.4/Progs/RNAcofold";
    $self->{RNAFILE} = undef;
    $self->{FOLDS}   = [];
    $self->{TMPFILE} = undef;
    $self->{RAWOUTPUT} = undef;

    bless $self;
    return $self;
   
}


sub getRawOutput {
  my ($self) = @_;
  return $self->{RAWOUTPUT};
}

sub setRNAFile {
    my ($self, $s) = @_;
    $self->{RNAFILE} = $s;
}


sub setCofoldSequences {
  my ($self, $seq, $mirna) = @_;

  $self->setRNASeq("$seq\&$mirna");
  
}


sub setRNASeq {
    my ($self, $s) = @_;
    
    my $tmpfile = Sets::getTempFile("/tmp/tota");
    open OUTFOLD, ">$tmpfile";
    print OUTFOLD ">rna\n$s\n";
    close OUTFOLD;
    
    $self->{RNAFILE} = $tmpfile;
    $self->{TMPFILE} = $tmpfile;
}


sub fold {
    my ($self) = @_;
    
    $self->{FOLDS} = [];
    
    my $out = `$self->{RNAFOLD} < $self->{RNAFILE}`;

    $self->{RAWOUTPUT} = $out;
    
    my @a = split /\n/, $out;

    my $s = undef;
    my $f = undef;
    my $e = undef;
    foreach my $l (@a) {
	
	next if ($l =~ /^\>/);
	
	if ($l =~ /[UACG]+/) {
	    $s = $l;
	} 

	if ($l =~ /^([\.\(\)\&]+)\ \(\ *([\-\d\.]+)\)/) {
	    $e = $2;
	    $f = $1;
	    
	    my %h = (ENERGY => $e, FOLD => $f, SEQUENCE => $s);
	
	    push @{ $self->{FOLDS} }, \%h;
	}
	
    
    

    }

    if (defined($self->{TMPFILE}) && -e $self->{TMPFILE}) {
      unlink $self->{TMPFILE};
    }

    return $self->{FOLDS};
}


1;
