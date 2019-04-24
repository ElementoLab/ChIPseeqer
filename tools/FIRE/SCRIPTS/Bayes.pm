package Bayes;

use strict ;
use Data::Dumper ;
use Sets;

sub new {
  my ($class) = @_;
  my ($self)  = {};

  $self->{FEATURE_NAMES}  = undef;
  $self->{FEATURE_MATRIX} = undef;
  $self->{CLASS_VECTOR}   = undef;
  $self->{OBJECT_NAMES}   = undef;
  $self->{PROBS}          = undef; # in  of cluster
  $self->{PROBS_C}        = undef;
  
  bless $self;
  return $self;

}


sub setFeatureNames {
  my ($self, $a_ref) = @_;
  $self->{FEATURE_NAMES} = $a_ref;
}

sub setFeatureMatrix {
  my ($self, $a_ref) = @_;
  $self->{FEATURE_MATRIX} = $a_ref;
}

sub setClassVector {
  my ($self, $a_ref) = @_;
  $self->{CLASS_VECTOR} = $a_ref;  
}

sub setObjectNames {
  my ($self, $a_ref) = @_;
  $self->{OBJECT_NAMES} = $a_ref;
}

sub train {
  my ($self) = @_;

  my $a_ref = Sets::transpose( $self->{FEATURE_MATRIX} );

  my $n_f = @$a_ref;
  
  my @counts_C = ();
  foreach my $r (@{$self->{CLASS_VECTOR}}) {
    $counts_C[ $r ] ++;
  }	

  $self->{PROBS_C} = log($counts_C[1])  - log ($counts_C[0]);

  my $j = 0;
  foreach my $f (@$a_ref) {
    #
    # calculate P(F=1|C=1) = P(F=1,C=1) / P(C=1)
    #           P(F=0|C=1)
    #           P(F=1|C=0)
    #           P(F=0|C=0)
    #
    
    my @counts = ();    
    my $maxf   = 0;
    for (my $i=0; $i<@$f; $i++) {
      $counts[$f->[$i]][ $self->{CLASS_VECTOR}->[$i]] ++;
      if ($f->[$i] > $maxf) {
	$maxf = $f->[$i];
      }
    }

    #print "Feature $j\n";
    # log P(F=1|C=1) - log P(F=1|C=0)
    for (my $k=0; $k<=$maxf; $k++) {
      $self->{PROBS}->[$j][$k] = 
	log ($counts[$k][1] / $counts_C[1] ) -
	  log ($counts[$k][0] / $counts_C[0] );
      #print " F=$k, log [ P(F=$k|C=1) / P(F=$k|C=0) ] = " . $self->{PROBS}->[$j][$k] . "\n";
    }
    

    $j++;
  }

}


sub classify {
  my ($self, $a_ref_F) = @_;

  print "Prior: $self->{PROBS_C}\n";
  
  my @classes = ();
  for (my $i=0; $i<@$a_ref_F; $i++) {
    my $lp = $self->{PROBS_C};
    for (my $j=0; $j<@{$a_ref_F->[$i]}; $j++) {
      $lp += $self->{PROBS}->[$j]->[ $a_ref_F->[$i]->[$j]  ];
    }
    if ($lp > 0) {
      $classes[$i] = 1;
    } else {
      $classes[$i] = 0;
    }
  }

  return \@classes;

}


1 ;
