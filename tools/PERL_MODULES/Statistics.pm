package Statistics;

#
# source: http://en.wikipedia.org/wiki/Student's_t-test
#
sub t_statistics {
  my ($a_ref1, $a_ref2) = @_;

  my $n1 = @$a_ref1;
  my $n2 = @$a_ref2;
  my $a1 = &average($a_ref1);
  my $a2 = &average($a_ref2);
  my $s1 = &stddev ($a_ref1);
  my $s2 = &stddev ($a_ref2);

  my $s  = sqrt( (1/$n1 + 1/$n2) * ( ($n1 - 1) * $s1 * $s1 + ($n2 - 1) * $s2 * $s2 ) / ($n1 + $n2  - 2) );
  
  my $t  = ($a1 - $a2) / $s;

  return $t;

}






sub average {
    my ($a_ref) = @_;
    
    my $m = scalar(@$a_ref);

    # calculate the sums
    my $sum   = 0.0;
    
    for (my $i=0; $i<$m; $i++) {
	$sum   += $a_ref->[$i];
	
    }
    return $sum / $m;
    
    
}




#
#  calculate the standard deviation of a series of values
#
sub stddev {
    my ($a_ref) = @_;
    
    my $m = scalar(@$a_ref);

    return -1 if ($m <= 1);

    # calculate the sums
    my $sum   = 0.0;
    my $sum_2 = 0.0;
    for (my $i=0; $i<$m; $i++) {
	$sum   += $a_ref->[$i];
	$sum_2 += $a_ref->[$i] * $a_ref->[$i];
    }
   
    # calculate the standard deviation
    my $std = sqrt( ($sum_2 - $sum * $sum / $m ) / ( $m - 1 )); 
  
    return $std;

}





1;
