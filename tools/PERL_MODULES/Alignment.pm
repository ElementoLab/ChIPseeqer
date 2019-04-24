package Alignment;
use strict;


#usage perl ~/PERL_MODULES/TEST/test_Alignment.pm

sub new {
  my ($class) = @_;

  my $self = {};
  $self->{P} =  -1; #-1;
  $self->{R} =   1; # 3;
  $self->{G} =  -1; #-1;
  $self->{O} =  -2; #-1;
  
  $self->{VERBOSE} = 0;
  $self->{COMP} = 'DNA';
  bless($self);
  return $self;
}


sub setVerbose {
  my ($self, $n) = @_;
  $self->{VERBOSE} = $n;
}

sub setComp {
  my ($self, $n) = @_;
  $self->{COMP} = $n;
}

sub rnamatch {
  my ($self, $n1, $n2) = @_;
  
  if ((($n1 eq 'A') && ($n2 eq 'T')) || (($n1 eq 'T') && ($n2 eq 'A'))) {
    return $self->{R};
  } elsif ((($n1 eq 'G') && ($n2 eq 'C')) || (($n1 eq 'C') && ($n2 eq 'G'))) {
    return $self->{R};
  } elsif ((($n1 eq 'G') && ($n2 eq 'T')) || (($n1 eq 'T') && ($n2 eq 'G'))) {
    return $self->{R};
  } else {
    return $self->{P};
  } 
  
  
  
}

sub dnamatch {
  my ($self, $n1, $n2) = @_;

  if ($n1 eq $n2) {
    return $self->{R};
  }  else {
    return $self->{P};
  }
  
}

sub score {
  my ($self, $n1, $n2) = @_;
  if ($self->{COMP} eq 'DNA') {
    return $self->dnamatch($n1, $n2);
  } elsif ($self->{COMP} eq 'RNA') {
    return $self->rnamatch($n1, $n2);    
  }
}
  
sub nw {
  my ($self, $ss1, $ss2) = @_;

  #print "$ss1\n";
  #print "$ss2\n";

  my $s1 = "N$ss1";
  my $s2 = "N$ss2";

  my @a1 = split //, $s1;
  my @a2 = split //, $s2;

  my @A  = ();  # scores
  my @L  = ();  # links

  # vertical
  for (my $i=0; $i<@a1; $i++) {
    $A[$i][0] = 0;
  }
 
  # horizontal
  for (my $j=0; $j<@a2; $j++) {
    $A[0][$j] = 0;
  }

  $L[0][0]->[0] = [ 0, -1, -1 ];
  
  for (my $i=0; $i<@a1; $i++) {
    
    for (my $j=0; $j<@a2; $j++) {

      next if (($i == 0) && ($j == 0));
      
      # printf "ENTRY i=$i, j=$j\n";

      # determine penalty or reward
      my $sc = $self->score($a1[$i], $a2[$j]);
      # print "$sc ($a1[$i] $a2[$j])\n";

      # now determine where the max is
      my @b = ();

      if ( ($i>0) && ($j>0) ) {
	my @a_tmp = ($A[$i-1][$j-1] + $sc, $i-1, $j-1, 0);
	push @b, \@a_tmp;
      }
      
      for (my $ii = $i-1; $ii>=0; $ii--) {
	my @a_tmp = ($A[$ii][$j] + ($i-$ii) * $self->{G} + $self->{O}, 
		     $ii, $j); 
	push @b, \@a_tmp;
      }
      
      for (my $jj = $j-1; $jj>=0; $jj--) {
	my @a_tmp = ( $A[$i][$jj] + ($j-$jj) * $self->{G} + $self->{O}, 
		      $i, $jj); 
	push @b, \@a_tmp;
      }
  
      @b = sort { $b->[0] <=> $a->[0] } @b; 
      
      #foreach my $s (@b) {
      #	print " [ $s->[0] ($s->[1], $s->[2]) ]";
      #}
      #print "\n";
      
      my $max = shift @b;

      $A[$i][$j] = $max->[0];


      push @{ $L[$i][$j] }, $max;
      
      my $cnt = 0;
      while ($b[$cnt] && ($b[$cnt]->[0] == $max->[0])) {
      	push @{ $L[$i][$j] }, $max;
      	$cnt ++;
      }
      
      

    }
    
  }
  
  if ($self->{VERBOSE} == 1) {
  
    print "    ";
    for (my $j=0; $j<@a2; $j++) {
      print sprintf("%4d", $j);
    }
    print "\n";
    
    for (my $i=0; $i<@a1; $i++) { 
      print sprintf("%4d", $i);
      for (my $j=0; $j<@a2; $j++) {
	print sprintf("%4d", $A[$i][$j]);
      }
      print "\n";
    }
    print "\n";
    
    for (my $i=0; $i<@a1; $i++) {    
      for (my $j=0; $j<@a2; $j++) {
	print sprintf(" [%2d %2d]", $L[$i][$j]->[0]->[1], $L[$i][$j]->[0]->[2]);
      }    
      print "\n";
    }
  }

  my $i  = @a1-1;
  my $j  = @a2-1;
  my @ai = ();
  my @aj = ();

  #print "$i\t$j\n";

  while ( 1 ) {

    my $l    = $L[$i][$j]->[0];

    my $myi  = $l->[1];
    my $myj  = $l->[2];

    last if (($myi==-1) && ($myj==-1));

    my $li   = $i - $myi;
    my $lj   = $j - $myj;

    #print "L $li\t$lj\n";

    if ($li == 0) {
      
      for (my $k=$j; $k>$myj; $k--) {
	push @ai, '-';
	push @aj, $a2[$k];
      }

    } elsif ($lj == 0) {

      for (my $k=$i; $k>$myi; $k--) {
	push @ai, $a1[$k];
	push @aj, '-';
      }

    } else {
      push @ai, $a1[$i];
      push @aj, $a2[$j];
    }
    
    $i = $myi;
    $j = $myj;

    #print "$i\t$j\n";
    
    #print "CUR " . join("", reverse(@ai)) . "\n";
    #print "CUR " . join("", reverse(@aj)) . "\n";
 
  }

  my @a =  ( join("", reverse(@ai)), join("", reverse(@aj)) );

  return \@a;
}




  
sub sw {
  my ($self, $ss1, $ss2) = @_;

  #print "$ss1\n";
  #print "$ss2\n";

  my $s1 = "N$ss1";
  my $s2 = "N$ss2";

  my @a1 = split //, $s1;
  my @a2 = split //, $s2;

  my @A  = ();  # scores
  my @L  = ();  # links

  my $be = ();  # keep track of best score and coords

  # vertical
  for (my $i=0; $i<@a1; $i++) {
    $A[$i][0] = 0;
  }
 
  # horizontal
  for (my $j=0; $j<@a2; $j++) {
    $A[0][$j] = 0;
  }

  $L[0][0]->[0] = [ 0, -1, -1 ];
  
  for (my $i=0; $i<@a1; $i++) {
    
    for (my $j=0; $j<@a2; $j++) {

      next if (($i == 0) && ($j == 0));
      
      # printf "ENTRY i=$i, j=$j\n";

      # determine penalty or reward
      my $sc = $self->score($a1[$i], $a2[$j]);
      # print "$sc ($a1[$i] $a2[$j])\n";

      # now determine where the max is
      my @b = ();

      if ( ($i>0) && ($j>0) ) {
	my @a_tmp = ($A[$i-1][$j-1] + $sc, $i-1, $j-1, 0);
	push @b, \@a_tmp;
      }
      
      for (my $ii = $i-1; $ii>=0; $ii--) {
	my @a_tmp = ($A[$ii][$j] + ($i-$ii) * $self->{G} + $self->{O}, 
		     $ii, $j); 
	push @b, \@a_tmp;
      }
      
      for (my $jj = $j-1; $jj>=0; $jj--) {
	my @a_tmp = ( $A[$i][$jj] + ($j-$jj) * $self->{G} + $self->{O}, 
		      $i, $jj); 
	push @b, \@a_tmp;
      }

      # smith-waterman 
      push @b, [0, -1, -1];
  
      @b = sort { $b->[0] <=> $a->[0] } @b; 
      
      #foreach my $s (@b) {
      #	print " [ $s->[0] ($s->[1], $s->[2]) ]";
      #}
      #print "\n";
      
      # maximum entry
      my $max = shift @b;


      # store maximum
      $A[$i][$j] = $max->[0];

      if ($max->[0] > $be->[0]) {
	my @nbe = ($max->[0], $i, $j);
	print "Max $max->[0] at $i, $j\n";
	$be = \@nbe;
      }

      
      # store link to max
      push @{ $L[$i][$j] }, $max;
      
      # and all other link to max
      my $cnt = 0;
      while ($b[$cnt] && ($b[$cnt]->[0] == $max->[0])) {
      	push @{ $L[$i][$j] }, $max;
      	$cnt ++;
      }
      
      

    }
    
  }
  
  if ($self->{VERBOSE} == 1) {
  
    print "    ";
    for (my $j=0; $j<@a2; $j++) {
      print sprintf("%4d", $j);
    }
    print "\n";
    
    for (my $i=0; $i<@a1; $i++) { 
      print sprintf("%4d", $i);
      for (my $j=0; $j<@a2; $j++) {
	print sprintf("%4d", $A[$i][$j]);
      }
      print "\n";
    }
    print "\n";
    
    for (my $i=0; $i<@a1; $i++) {    
      for (my $j=0; $j<@a2; $j++) {
	print sprintf(" [%2d %2d]", $L[$i][$j]->[0]->[1], $L[$i][$j]->[0]->[2]);
      }    
      print "\n";
    }
  }

  # where we start (max score) 
  my $i  = $be->[1]; #@a1-1;
  my $j  = $be->[2]; #@a2-1;
  my @ai = ();
  my @aj = ();

  print "START: $i\t$j\n";

  my @st = ($i, $j);
  my @en = ();

  while ($A[$i][$j] > 0) {

    my $l    = $L[$i][$j]->[0];

    

    my $myi  = $l->[1];
    my $myj  = $l->[2];

    if ($self->{VERBOSE} == 1) {
      printf("Go to $myi, $myj\n");
    }	
    last if (($myi==-1) && ($myj==-1));

    my $li   = $i - $myi;
    my $lj   = $j - $myj;

    #print "L $li\t$lj\n";

    if ($li == 0) {
      
      for (my $k=$j; $k>$myj; $k--) {
	push @ai, '-';
	push @aj, $a2[$k];
      }

    } elsif ($lj == 0) {

      for (my $k=$i; $k>$myi; $k--) {
	push @ai, $a1[$k];
	push @aj, '-';
      }

    } else {
      push @ai, $a1[$i];
      push @aj, $a2[$j];
    }
    
    $i = $myi;
    $j = $myj;

    @en = ($i, $j);
    
    #print "CUR " . join("", reverse(@ai)) . "\n";
    #print "CUR " . join("", reverse(@aj)) . "\n";
 
  }

  my @a =  ( [join("", reverse(@ai)), $en[0], $st[0]], [join("", reverse(@aj)),$en[1], $st[1]] );

  return \@a;
}





1;
