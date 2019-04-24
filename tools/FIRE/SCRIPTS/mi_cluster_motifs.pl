use lib "$ENV{FIREDIR}/SCRIPTS";

use AggloClust;
use Table;
use Sets;
use strict;
# use Data::Dumper;

my $nbclusters   = Sets::get_parameter(\@ARGV, "-nbclusters");
my $matrixfile   = Sets::get_parameter(\@ARGV, "-matrixfile");
my $clusterfile  = Sets::get_parameter(\@ARGV, "-clusterfile");

# read in matrix file
my $ta = Table->new;
$ta->loadFile($matrixfile);
my $a_ref_mat = $ta->getArray();
shift @$a_ref_mat;
 
my $ac = AggloClust->new;

my @dist = ();
my $n    = @$a_ref_mat;
for (my $i=0; $i<$n-1; $i++) {
  $dist[$i][$i] = 0;
  for (my $j=$i+1; $j<$n; $j++) {
    my @a1 = @{ $a_ref_mat->[$i] }; shift @a1; 
  
#    foreach my $r (@a1) {
#      if ($r > Sets::log10(0.05/@a1)) {
#	$r =  1;
#      } elsif ($r > Sets::log10(0.05/@a1)) {
#	$r =  0;
#      } else {
#	$r =  0;
#      }
#    }

    my @a2 = @{ $a_ref_mat->[$j] }; shift @a2; 

#    foreach my $r (@a2) {
#      if ($r > Sets::log10(0.05/@a2)) {
#	$r =  1;
#      } elsif ($r > Sets::log10(0.05/@a2)) {
#	$r =  0;
#      } else {
#	$r =  0;
#      }
#    }


    $dist[$i][$j] = 1 - Sets::pearson(\@a1, \@a2);      
    $dist[$j][$i] = $dist[$i][$j]; 
  }
}

$ac->setDistanceMatrix(\@dist);
#$ac->setMaxNbClusters($nbclusters);

my $a_ref_c = $ac->agglomerate_using_max_linkage();

my $lastn = @$a_ref_c;

print "LASTN $lastn\n";

for (my $i=0; $i<$lastn; $i++) {
  print "$i => [ " . $a_ref_c->[$i]->[0] . " " . $a_ref_c->[$i]->[1] . " ]\n"; 
}

#exit;

#print Dumper($a_ref_c);

&visit($lastn-1, $a_ref_c);
  



sub visit {
  my ($node, $tree) = @_;

  print "Visiting $node\n";

  my ($l, $r) = @{ $tree->[$node] };

  #print "Children $l and $r\n";


  if (!defined($l) && !defined($r)) {
    
    print "Reached leaf $node\n";

    return;
  }

  
  
  &visit($l, $tree);
  
  
  &visit($r, $tree);
}


exit;

open OUT, ">$clusterfile";

my $cnt = 0;
foreach my $c (@$a_ref_c) {
  print join(" ", @$c); print "\n";
  foreach my $i (@$c) {
    print OUT "$a_ref_mat->[$i]->[0]\t$cnt\n";
  }
  $cnt ++;
}
  
close OUT;
