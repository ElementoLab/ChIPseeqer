package Graph;

sub new {
  
  my $self = {};

  # adjacency list
  $self->{ADJ}         = {};  

  bless($self);
  return $self;
  
}

#
# load graph (list of edges)
#
sub loadGraph {
  my ($self, $file) = @_;

  open IN, $file;
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;

    push @{ $self->{ADJ}->{ $a[0] } }, $a[1];
    push @{ $self->{ADJ}->{ $a[1] } }, $a[0];

  }
  close IN;

}



sub getCenteredGraph {
  my ($self, $g, $level) = @_;

  # this is a DFS, with restricted levels

  # init colors
  my %COL = ();
  foreach my $n (keys(%{ $self->{ADJ} })) {
    $COL{$n} = 0;
  }

  my %COL_EDGES = ();

  my @a_edges = ();
  $self->_DFS_edges($g, $g, \%COL, \%COL_EDGES, \@a_edges, 0, $level);

  return \@a_edges;

}



#
#  returns connected components out of an adjacency list H{ node } => node
#
sub getConnectedComponents {
  my ($self) = @_;
  
  my @COMP = ();

  # init colors
  my %COL = ();
  foreach my $n (keys(%{ $self->{ADJ} })) {
    $COL{$n} = 0;
  }
  
  # DFS
  foreach my $n (keys(%{ $self->{ADJ} })) {
    if ($COL{$n} == 0) {
      my @a_nodes = ();
      _DFS($n, \%COL, \@a_nodes);
      push @COMP, \@a_nodes;
    }
  }

  return \@COMP;
} 
	

sub _DFS {
    my ($self, $node, $h_ref_col, $a_ref_nodes, $level, $levelmax) = @_;

    $h_ref_col->{ $node } = 1;

    if (defined($level) && defined($levelmax)&& ($level == $levelmax)) {
      return;
    }

    foreach my $next (@{ $self->{ADJ}->{ $node } }) {
	next if (!defined($self->{ADJ}->{$next}));
	if ($h_ref_col->{ $next } == 0) {
            _DFS($next, $h_ref_col, $a_ref_nodes);
        }
    }
    $h_ref_col->{ $node } = 2;
}


#
# DFS that returns all edges
#
sub _DFS_edges {

  my ($self, $node, $prevnode, $h_ref_col, $h_ref_col_edges, $a_ref_nodes, $level, $levelmax) = @_;

  
  
  #print "Entering node $node\n";
  
  # adding edge
  if (($prevnode ne $node) && (!defined($h_ref_col_edges->{$prevnode}->{$node}))) {

    push @$a_ref_nodes, [ $prevnode, $node ];

    $h_ref_col_edges->{$prevnode}->{$node} = 1;
    $h_ref_col_edges->{$node}->{$prevnode} = 1;

  }


  if ($h_ref_col->{ $node } > 0) {
    # already visited
    return;
  }

  $h_ref_col->{ $node } = 1;


  if (defined($level) && defined($levelmax)&& ($level == $levelmax)) {
    return;
  }

  foreach my $next (@{ $self->{ADJ}->{ $node } }) {
    next if (!defined($self->{ADJ}->{$next}));

    $self->_DFS_edges($next, $node, $h_ref_col, $h_ref_col_edges, $a_ref_nodes, $level+1, $levelmax);
    
  }
  



}





1;
