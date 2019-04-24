package Markov;

use strict;

sub new {
    my $self  = {};

    $self->{ ALP} = ['A', 'C', 'T', 'G'];
    $self->{NALP} = ['A', 'C', 'T', 'G', 'N'];

    $self->{USEN} = 0;
    
    # single nucleotide counts
    $self->{COUNTS1}  = undef;
    $self->{SIZE1}    = undef;

    # dinucleotides freq
    $self->{COUNTS2}  = undef;
    $self->{SIZE2}    = undef;
    
    foreach my $l (@{ $self->{NALP} }) {
	$self->{COUNTS1}->{$l} = 0.01;
    }
    
    foreach my $l1 (@{ $self->{NALP} }) {
	foreach my $l2 (@{ $self->{NALP} }) {
	    $self->{COUNTS2}->{$l1}->{$l2} = 0.01;
	}
    }

    $self->{SIZE1}   = 0;
    $self->{SIZE2}   = 0;
    $self->{VERBOSE} = 0;
    
    srand();

    bless($self);           
    return $self;
}


sub useN {
   my ($self, $i) = @_;
   $self->{ USEN } = $i;
}



sub calcFrequenciesFromSeq {
    my ($self, $s_str) = @_;
    
    my $h_ref_tmp = $self->countNucleotides($s_str);	

    $self->_addCount($h_ref_tmp); #length($s_str));
    
    #$m1->printCounts;


}




# for a specific sequence
sub countNucleotides {
    my ($self, $s_seq) = @_;

    my %h_cntbkg = ('A'     => 0,
		    'C'     => 0,
		    'T'     => 0,
		    'G'     => 0,
		    'N'     => 0,
		    'LEN'   => 0,
		    'DILEN' => 0);
    
    #
    # initialise the 2-nt count
    #
    foreach my $l1 (@{ $self->{NALP} }) {
      foreach my $l2 ( @{ $self->{NALP} } ) {
	$h_cntbkg{$l1 . $l2} = 0;
      }
    }
    
    
    # analyze the sequence
    my @a_seq    = split //, $s_seq;
    for (my $i=0; $i < length($s_seq); $i++) {
      
      
      if ($self->{USEN} == 1) {
	next if ($a_seq[$i] !~ /[NATCG]/);
      } else {
	next if ($a_seq[$i] !~ /[ATCG]/);
      }

      # NT
      $h_cntbkg{$a_seq[$i]} ++;
      
      # INCREASE LENGTH
      $h_cntbkg{"LEN"} ++;

      #
      # no need to count dint if already at the end	
      #
      last if ($i == length($s_seq) - 1);
	
      if ($self->{USEN} == 1) {
	next if ($a_seq[$i+1] !~ /[NATCG]/);
      } else {
	next if ($a_seq[$i+1] !~ /[ATCG]/);
      }
      
      # DINT
      $h_cntbkg{$a_seq[$i] . $a_seq[$i+1]} ++;
      
      # INCREASE LENGTH
      $h_cntbkg{"DILEN"} ++;

    }

    return \%h_cntbkg;
}


#
#  add to the current counts
#
sub _addCount {
    my ($self, $h_ref) = @_;
    
    foreach my $l ( @{ $self->{NALP} } ) {
      $self->{COUNTS1}->{$l} += $h_ref->{$l};
    }
    
    foreach my $l1 ( @{ $self->{NALP} } ) {
      foreach my $l2 ( @{ $self->{NALP} } ) {
	$self->{COUNTS2}->{$l1}->{$l2} += $h_ref->{$l1 . $l2};
      }
    }
    
    $self->{SIZE1} += $h_ref->{  LEN};
    $self->{SIZE2} += $h_ref->{DILEN};
    
}



# generate a 1-order sequence of length L
sub generate1rstOrder {
    my ($self, $i_len) = @_;

    my $l =  $self->_generateFirstLetter();

    my $s_seq = $l;

    for (my $i=1; $i<$i_len; $i++) {
	$l      = $self->_generateNextLetter($l);
	$s_seq .= $l;
    }
    
    die "len(s_seq) != $i_len" if (length($s_seq) != $i_len);

    return $s_seq;
}

#
# generate a first letter from a Markov model
#
sub _generateFirstLetter {
  my ($self) = @_;

  # create an array of cumulative probabilities
  my @a_cum = ();
  my $i     = 1;

  my @a     = undef;  
  if ($self->{USEN} == 1) {
    @a = @{ $self->{NALP} };
  } else {
    @a = @{ $self->{ ALP} };
  }
  my $m     = scalar( @a );
  
  $a_cum[0] = 0.0;
  foreach my $l1 (@a) {
    $a_cum[$i] = $a_cum[$i - 1] + $self->{COUNTS1}->{$l1} / $self->{SIZE1};
    $i++;
  }
  
  my $d = rand;
  my $n = undef;
  
  for ($i=1; $i<=$m; $i++) {
    if (($d > $a_cum[$i-1]) && ($d <= $a_cum[$i])) {
      $n = $a[$i-1];
      last;
    }
  }
  
  
  return (defined($n)?$n:'N');
    
}


sub _generateNextLetter {

    my ($self, $l) = @_;

    my @a_cum = ();
    my $i     = 1;

    my @a     = undef;  
    if ($self->{USEN} == 1) {
      @a = @{ $self->{NALP} };
    } else {
      @a = @{ $self->{ ALP} };
    }
    my $m     = scalar( @a );

    #
    #  create a cumulative table for the current 
    #
    $a_cum[0] = 0.0;
    foreach my $l1 (@a) {
	$a_cum[$i] = $a_cum[$i - 1] + $self->{COUNTS2}->{$l}->{$l1} / $self->{SIZE2};
	$i++;
    }

    # normally, index $i-1 should be equal to one, so normalize all values
    my $t = 1.0 / $a_cum[$i - 1];
    
    for ($i=0; $i<=$m; $i++) {
	$a_cum[$i] = $a_cum[$i] * $t;
    }
        
    # draw a float between 0 and 1
    my $d = rand;
    
    #print "Drawn $d ..\n";
    my $n = undef;
    
    for ($i=1; $i<=$m; $i++) {
      if (($d > $a_cum[$i-1]) && ($d <= $a_cum[$i])) {
	$n = $a[$i-1];
	last;
      }
    }
    

    return (defined($n)?$n:'N');
    
    
}






sub printCounts {
    my ($self) = @_;
	

    #print "SIZE=$self->{SIZE}\n";
    my $i1 = 0;
    my %h_tmp = %{$self->{COUNTS1}};
    reset %h_tmp;    
    while (my ($l, $n) = each(%h_tmp)) {
	print sprintf("%s\t%d\t%5.4f\n", $l, $n, $n/$self->{SIZE});
	$i1 += $n;
    }

     my $i2 = 0;
    my %h_tmp1 = %{$self->{COUNTS2}};
    reset %h_tmp1;    
    while (my ($l1, $n1) = each(%h_tmp1)) {
	
	my %h_tmp2 = %$n1;	    
	reset %h_tmp2;    
	 while (my ($l2, $n2) = each(%h_tmp2)) {
	     print sprintf("%s\t%d\t%5.4f\n", $l1 . $l2, $n2, $n2/$self->{SIZE});
	     	$i2 += $n2;
	 }
    }

    #if ($i1-1 != $i2 ) {
#	die "$i1-1 != $i2 .. !\n";
    #}


    
}




# 
sub writeFrequencies {
    
    my ($self, $s_file) = @_;

    open OUT, ">$s_file";


    my %h_tmp = %{$self->{COUNTS1}};
    reset %h_tmp;    
    while (my ($l, $n) = each(%h_tmp)) {
	print OUT sprintf("%s\t%f\n", $l, $n/$self->{SIZE});
    }

    my %h_tmp1 = %{$self->{COUNTS2}};
    reset %h_tmp1;    
    while (my ($l1, $n1) = each(%h_tmp1)) {
	my %h_tmp2 = %$n1;	    
	reset %h_tmp2;    
	 while (my ($l2, $n2) = each(%h_tmp2)) {
	     print OUT sprintf("%s\t%f\n", $l1 . $l2, $n2/$self->{SIZE});
	 }
    }

    close OUT;

}

# generate a 0-order Markov sequence of length L
sub generateOthOrder {
    my ($self, $l) = @_;

    


}


sub getFrequencies {
    my ($self) = @_;

    my %h = ();
    reset(%{$self->{COUNTS1}});
    while (my ($k, $v) = each(%{$self->{COUNTS1}})) {
	
	#print "h($k) = " . $v / $self->{SIZE} . "\n";
	$h{$k} = $v / $self->{SIZE};
	
    }

    
    return \%h;
}


# get the dinucleotide frequencies
sub getDiFrequencies {
    
    my ($self) = @_;

    my %h = ();
    reset(%{$self->{COUNTS2}});
    while (my ($k1, $v1) = each(%{$self->{COUNTS2}})) {
	
	reset(%$v1);
	while (my ($k2, $v2) = each(%$v1)) {


	    #print "h($k1$k2) = " . $v2 / $self->{SIZE} . "\n";
	    $h{$k1 . $k2} = $v2 / $self->{SIZE};
	    
	}
	
    }

    return \%h;
    
}



# read from file
sub readFrequencies {
    
    my ($self, $s_file) = @_;

    open IN, $s_file or die "Could not open file $s_file ...\n";
    while (my $s_line = <IN>) {
	
	chomp $s_line;
	
	my @a_tmp = split /\t/, $s_line;
	
	if (length($a_tmp[0]) == 1) {
	    $self->{COUNTS1}->{$a_tmp[0]} = $a_tmp[1];
	    
	} elsif (length($a_tmp[0]) == 2) {
	    my ($l1, $l2) = split //, $a_tmp[0];
	    $self->{COUNTS2}->{$l1}->{$l2} = $a_tmp[1];
	}
	
    }

    $self->{SIZE} = 1.0;

}

1;

