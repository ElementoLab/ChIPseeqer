package Log;

sub new {
    my ($class) = @_;
    my ($self) = {};
    $self->{FILE}  = 'log.txt';

    open LOG, ">log.txt";

    bless $self;
    return $self;

    
}


sub log {
  my ($self, $txt) = @_;
    
  print LOG $txt; 
    

}

sub DESTROY {
  my ($self) = @_;

  close LOG;
}


1;
