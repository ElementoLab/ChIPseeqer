BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

#
# read in windows
#

my %POS = ();

open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @{ $POS{ $a[0] } }, $a[ 2 ] if (!Sets::in_array($a[2], @{ $POS{ $a[0] } }));
}
close IN;


print "loaded positions\n";

foreach my $k (keys(%POS)) {
  @{ $POS{ $k } } = sort { $a <=> $b } @{ $POS{ $k } };
  &output_windows(\@{ $POS{ $k } }, 500, 3, $k);
}




sub output_windows{
  
  my ($a_ref_positions, $w, $k, $chr) = @_; 

  my $n = scalar(@$a_ref_positions);
  
  my $global_k;
  my $end_j;
  my $max_local_k;
  my $local_k;
  my $sum_pos;
  my $max_local_k_avg_pos;


  for (my $i=0; $i<$n-1; $i++) {
    
    $global_k = 1;


    $end_j = $i;  
    for ($j=$i+1; $j<$n; $j++) {      
      if ( ($a_ref_positions->[ $j ] - $a_ref_positions->[ $j-1 ]) <= $w) {
	$end_j = $j;
	$global_k ++;
      } else {
	last;
      }
    }

    
    if ($global_k >= $k) {
      
      $max_local_k = -1;

      for ($l=$i; $l<=$end_j; $l++) {
	
	$sum_pos = $a_ref_positions->[$l];
	$local_k = 1;

	for ($m=$l+1; $m<=$end_j; $m++) {

	  if ( ($a_ref_positions->[$m] - $a_ref_positions->[$l]) <= $w) {
	    $sum_pos += $a_ref_positions->[$m];
	    $local_k ++;
	  } else {
	    last;
	  }

	}
	
	if (($local_k >= $k) && ($local_k > $max_local_k)) {
	  $max_local_k         = $local_k;
	  $max_local_k_avg_pos = $sum_pos / $local_k;
	}
      }
      
      if ($max_local_k > 0) {
	printf("%s\t%d\t%d\n", $chr, $max_local_k_avg_pos, $max_local_k);
      }
      
    }
   
    $i = $end_j;
 
  }

}
		
