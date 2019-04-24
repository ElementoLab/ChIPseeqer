# count the number of clusters for which there is a 1-to-1 rekationship
BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Sets;
use Hypergeom;
use strict;


# read in blastclust clusters
open IN, $ARGV[0];
my @C1 = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l; 
  while ($a[ $#a ] eq "") {
    pop @a;
  }

  #
  #  make sure that @a has at least 5 species
  #
  
  my %HCNT = ();
  foreach my $s (@a) {
    my @o = split /\.1_/, $s;
    $HCNT { $o[1] } ++;  
  }
  
    
  my $cnt_sp_max = scalar(keys(%HCNT));

  next if ($cnt_sp_max < 5);
  

  push @C1, \@a;
}
close IN;


# read in GG file names
my $a_ref_f = Sets::getFiles($ARGV[1]);
my $nc2 = scalar(@$a_ref_f);
# read in GG cluster
my @C2 = ();
foreach my $f (@$a_ref_f) {
  
  my @C = ();
  open IN, $f;
  while (my $l = <IN>) {
    chomp $l;
    next if (($l =~ /^\%/) || ($l eq ""));
    my @a = split /\t/, $l; 
    next if ($a[1] == 0);
    my $g = $a[3] . "_" . $a[5];
    push @C, $g;
  }

  push @C2, \@C;
  
  close IN;
  
}


# now run the comparison ?

my $recip  = 0;
my $nrecip = 0;
my $p_best = -1;
my $i      = 0;
foreach my $c1 (@C1) {
  
  my $i2 = undef;
  my $p2 = undef;
  &find_index_of_strongest_overlapping_cluster($c1, \@C2, $ARGV[2], \$i2, \$p2);
  
  # print "closest to c1[ $i ] in C2 is $i2\n"; 

  my $i1 = undef;
  my $p1 = undef;
  &find_index_of_strongest_overlapping_cluster($C2[$i2], \@C1, $ARGV[2], \$i1, \$p1);

  if ($i == $i1) {
    $recip ++; if ($p1 > $p_best) { $p_best = $p1; }
    #print "C1[$i] => C2[$i2], log(p)=" . sprintf("%4.3e", $p1) . "\n";
  } else {
    $nrecip ++;
  }

  $i ++;
  
}

print "recip = $recip, nrecip = $nrecip, Nggs = $nc2, best p = $p_best\n";


sub find_index_of_strongest_overlapping_cluster {
  my ($a_ref_c1, $a_ref_C2, $N, $ref_i, $ref_p) = @_;

  my $n = scalar(@$a_ref_C2);

  my $i_best = 0;
  my $p_best = 10.0;
  for (my $i=0; $i<$n; $i++) {
    
    my $a_ref_ov = Sets::getOverlapSet($a_ref_c1, $a_ref_C2->[$i]);

    my $ov       = scalar(@$a_ref_ov);
    my $s1       = scalar(@$a_ref_c1);
    my $s2       = scalar(@{ $a_ref_C2->[$i] });
    
    my $p        = Hypergeom::cumhyper($ov, $s1, $s2, $N);
    
    if ($p < $p_best) {
      $p_best = $p;
      $i_best = $i;
    }
    
  }

  $$ref_i = $i_best;
  $$ref_p = $p_best;

  
}

