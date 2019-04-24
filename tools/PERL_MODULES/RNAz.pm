package RNAz;

sub new {
    my ($self)  = {};
    $self->{ISRNA} = undef;
 
    bless($self);
    return $self;
}

sub run {
  my ($self, $f) = @_;
  
  die "$f does not exist\n" if (! -e $f);
  
  my $todo = "RNAz < $f";
  my $out  = `$todo`;
  
  # parse the output file
  my @a = split /\n/, $out;

  foreach my $l (@a) {    
    if ($l =~ /Prediction\: (.+)$/) {
      if ($1 eq "no RNA") {
	$self->{ISRNA} = 0;
      } elsif ($1 eq "RNA") {
	$self->{ISRNA} = 1;
      } else {
	die "parsing problem\n";
      }
    }
  }
  
}

sub isRNA {
  my ($self) = @_;
  return $self->{ISRNA};
}


1;

