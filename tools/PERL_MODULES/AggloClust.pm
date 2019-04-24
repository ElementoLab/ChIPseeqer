package AggloClust;
use Sets;
use strict;

sub new {
  my ($class)         = @_;
  my ($self)          = {};
  
  $self->{DIST}       = [];
  $self->{MAXNBCLUST} = undef;
  $self->{MIND}       = undef;
  $self->{USECORR}    = 0;
  $self->{USEABS}     = 0;
  $self->{TREE}       = undef;
  $self->{ALGOCLUST}  = "max";
  
  bless $self;
  return $self;
}


sub setAlgoClust {
  my ($self, $m) = @_;
  $self->{ALGOCLUST} = $m;
}


sub setUseCorr {
  my ($self, $m) = @_;
  $self->{USECORR} = $m;
}

sub setUseAbs {
  my ($self, $m) = @_;
  $self->{USEABS} = $m;
}


sub setMaxNbClusters {
  my ($self, $m) = @_;
  $self->{MAXNBCLUST} = $m;
}


sub setDistanceMatrix {
  my ($self, $m) = @_;
  $self->{DIST} = $m;
}



sub setMinD {
  my ($self, $m) = @_;
  $self->{MIND} = $m;
}


sub getDFSOrder {
  my ($self) = @_;

  foreach  my $r (@{ $self->{TREE} }) {
    if ($r) {
      print join("\t", @$r) . "\n";
    } else {
      print "undef\n";
    }
  }

  my @order = ();
  $self->visit($self->{ROOTNODE}, $self->{TREE}, \@order);

  return \@order;
}


sub visit {
  my ($self, $node, $tree, $a_ref_order) = @_;

  # print "Visiting $node\n";

  my ($l, $r) = @{ $tree->[$node] };



  if (($l == -1) && ($r == -1)) {

    push @$a_ref_order, $node;
    
    # print "Reached leaf $node\n";

    return;
  }

  
  
  $self->visit($l, $tree, $a_ref_order);    
  $self->visit($r, $tree, $a_ref_order);
}


sub agglomerate_using_avg_linkage {

  my ($self) = @_;

  my $a_ref_d = $self->{DIST};
 
  if (!defined($self->{MAXNBCLUST}) && !defined($self->{MIND})) {
    die "Please set a stopping condition.\n";
  }

  my $n       = @$a_ref_d;
  
  #print "Got $n objects\n";
  
  my @CLASSES = ();
  my @LIST    = ();
  for (my $i=0; $i<$n; $i++) {
    my @a = ($i);
    push @CLASSES, \@a;
    push @LIST   , $i;
  }
 
  my @DIST = ();
  for (my $i=0; $i<$n; $i++) {
    for (my $j=0; $j<$n; $j++) {
      $DIST[$i][$j] = $a_ref_d->[$i]->[$j];
    }
  }

  my @NEWLIST = ();

  while (1) {
    my $a_ref_ijd = $self->find_min( \@DIST, \@LIST );
    my ($i, $j, $min) = @$a_ref_ijd;

    #print "best merge is $i, $j, $min\n";

    if ((defined($self->{MIND}) && ($min > $self->{MIND})) || (defined($self->{MAXNBCLUST}) && (scalar(@LIST) == $self->{MAXNBCLUST}))) {
      my @NEWCLASSES = ();
      foreach my $k (@LIST) {
	push @NEWCLASSES, $CLASSES[$k];
      }

      #
      # order the classes here ? first class is tightest ? (then use prior ??) P(X)
      #

      #print "order classes\n";

      my $best_myi = undef; my $best_d = 1000000;
      my $myi      = 0;
      foreach my $r (@NEWCLASSES) {
	if (scalar(@$r) > 1) {
	  my $d = $self->get_pairwise_avg_dist($r, $r) / (scalar(@$r) / $n);
	  #print scalar(@$r) . " - " . $d . "\n";
	  if ($d < $best_d) {
	    $best_myi = $myi;
	    $best_d   = $d;
	  }
	}
	$myi ++;
      }
      
      my $tmpa = $NEWCLASSES[0];
      $NEWCLASSES[0] = $NEWCLASSES[$best_myi];
      $NEWCLASSES[$best_myi] = $tmpa;

      return \@NEWCLASSES;
    }

    # merge i and j into i
    @{$CLASSES[$i]} = ( @{$CLASSES[$i]}, @{$CLASSES[$j]} );

    #print "Merge into " . join(" " , @{$CLASSES[$i]}) . "\n";

    # remove $j from list of classes
    @NEWLIST = ();
    foreach my $k (@LIST) {
      if ($k != $j) {
	push @NEWLIST, $k;
      }
    }
    @LIST = @NEWLIST;

    # update dist from new i to all other classes
    foreach my $k (@LIST) {
      next if ($k == $i); 
      $DIST[$i][$k] = $self->get_pairwise_avg_dist($CLASSES[$i], $CLASSES[$k]);
      #print "Updated dist[$i][$k] = $DIST[$i][$k]\n";
      $DIST[$k][$i] = $DIST[$i][$k];
    } 

    #<STDIN>;
  }

}




sub agglomerate_using_max_linkage {

  my ($self)   = @_;

  my $a_ref_d  = $self->{DIST};

  my $n        = @$a_ref_d;
  my @CLASSES  = ();
  my @LIST     = ();
  for (my $i=0; $i<$n; $i++) {
    my @a = ($i);
    push @CLASSES, \@a;
    push @LIST   , $i;
  }
 
  #
  #  copy distance matrix
  #
  my @DIST     = ();
  for (my $i=0; $i<$n; $i++) {
    for (my $j=0; $j<$n; $j++) {
      $DIST[$i][$j] = $a_ref_d->[$i]->[$j];
    }
  }
  
  my @NEWLIST = ();
  
  #
  # index of next entry in DIST
  #
  my $next_i = @DIST;

  # keep TREE structure
  my @NEXT = ();
  for (my $i=0; $i<$next_i; $i++) {
    $NEXT[$i] = [ -1, -1, -1 ];
  }

  while (@LIST > 1) {

    #
    # @LIST contains the list of indices in @DIST to visit
    #
    my $a_ref_ijd = undef;

    if ($self->{USECORR} == 0) {
      $a_ref_ijd = $self->find_min( \@DIST, \@LIST );
    } else {
      $a_ref_ijd = $self->find_max( \@DIST, \@LIST );
    }

    my ($i, $j, $min) = @$a_ref_ijd;

    #print "best merge is $i, $j, $min\n";

    # merge objects from classes i and j into next_i
    @{$CLASSES[$next_i]} = ( @{$CLASSES[$i]}, @{$CLASSES[$j]} );

    # add new class
    push @LIST, $next_i;

    $NEXT[$next_i] = [ $i, $j, $min ];

    # remove $i and $j from list of classes
    @NEWLIST = ();
    foreach my $k (@LIST) {
      if (($k != $i) && ($k != $j)) {
	push @NEWLIST, $k;
      }
    }
    @LIST = @NEWLIST;

    # update dist from next_i to all other classes
    foreach my $k (@LIST) {
      next if ($k == $next_i); 
      my $tmp = undef;

      if ($self->{ALGOCLUST} eq "max") {

	if ($self->{USECORR} == 0) {
	  $tmp = $self->get_pairwise_max_dist($CLASSES[$next_i], $CLASSES[$k]);
	} else {
	  $tmp = $self->get_pairwise_min_cor ($CLASSES[$next_i], $CLASSES[$k]);
	}

      } elsif ($self->{ALGOCLUST} eq "min") {
	
	if ($self->{USECORR} == 0) {
	  $tmp = $self->get_pairwise_min_cor ($CLASSES[$next_i], $CLASSES[$k]);
	} else {
	  $tmp = $self->get_pairwise_max_dist($CLASSES[$next_i], $CLASSES[$k]);
	}

      } elsif ($self->{ALGOCLUST} eq "avg") {
	
	# no need to chose between cor and non-corr
	$tmp = $self->get_pairwise_avg_dist($CLASSES[$next_i], $CLASSES[$k]);
	
      }
      
      $DIST[$next_i][$k] = $tmp;
      $DIST[$k][$next_i] = $DIST[$next_i][$k];
    }

    $next_i ++;

  }

  $self->{TREE}     = \@NEXT;
  $self->{ROOTNODE} = scalar(@NEXT) - 1;

  #return \@NEXT;

}

sub getTree {
  my ($self) = @_;

  return $self->{TREE};
  
}


#
#  calculate minimum pairwise correlation among two groups
#
sub get_pairwise_min_cor {
  my ($self, $a1, $a2) = @_;
  
  my $a_ref_d = $self->{DIST};

  my $n1 = @$a1;
  my $n2 = @$a2;
  
  my $sum = 10000000;

  for (my $i=0; $i<$n1; $i++) {
    for (my $j=0; $j<$n2; $j++) {

      if ($a_ref_d->[ $a1->[$i] ]->[ $a2->[$j] ] < $sum) {
	$sum = $a_ref_d->[ $a1->[$i] ]->[ $a2->[$j] ];
      }

    }
  }
  
  return $sum;
}



sub get_pairwise_max_dist {
  my ($self, $a1, $a2) = @_;
  
  my $a_ref_d = $self->{DIST};

  my $n1 = @$a1;
  my $n2 = @$a2;
  
  my $sum = 0;

  for (my $i=0; $i<$n1; $i++) {
    for (my $j=0; $j<$n2; $j++) {

      if ($a_ref_d->[ $a1->[$i] ]->[ $a2->[$j] ] > $sum) {
	$sum = $a_ref_d->[ $a1->[$i] ]->[ $a2->[$j] ];
      }

    }
  }
  
  return $sum;
}


sub get_pairwise_avg_dist {
  my ($self, $a1, $a2) = @_;
  
  my $a_ref_d = $self->{DIST};

  my $n1 = @$a1;
  my $n2 = @$a2;
  
  my $sum = 0;

  for (my $i=0; $i<$n1; $i++) {
    for (my $j=0; $j<$n2; $j++) {
      $sum += $a_ref_d->[ $a1->[$i] ]->[ $a2->[$j] ];
    }
  }
  
  $sum = $sum / ($n1 * $n2);
  
  return $sum;
}

sub find_min {
  my ($self, $a_ref_d, $a_ref_l) = @_;


  my $min     = 100000;
  my $min_i   = undef;
  my $min_j   = undef;

  foreach my $i (@$a_ref_l) { 
    foreach my $j (@$a_ref_l) { 
      next if ($i == $j);
      if ($a_ref_d->[$i]->[$j] < $min) {
	$min   = $a_ref_d->[$i]->[$j];
	$min_i = $i;
	$min_j = $j;
      }
    }
  }
  
  my $a = [$min_i, $min_j, $min];
  return $a;
}

#
#  finds best correlation in a correlation matrix
#
sub find_max {
  my ($self, $a_ref_d, $a_ref_l) = @_;


  my $max     = -100000;
  my $max_i   = undef;
  my $max_j   = undef;

  foreach my $i (@$a_ref_l) { 
    foreach my $j (@$a_ref_l) { 
      next if ($i == $j);
      my $t = $a_ref_d->[$i]->[$j];
      if ($self->{USEABS} == 1) {
	$t = abs($a_ref_d->[$i]->[$j]);
      }
      if ($t > $max) {
	$max   = $t;
	$max_i = $i;
	$max_j = $j;
      }
    }
  }
  
  my $a = [$max_i, $max_j, $max];
  return $a;
}


1;
