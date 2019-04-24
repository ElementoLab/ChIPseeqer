#
# this script is a very basic reimplementation of  
# Cancer Outlier Profile Analysis (COPA)
# from Tomlins et al, 2005.
# 
#   
#

if (@ARGV == 0) {
  die "Usage: perl COPA matrix.txt [0/1] (matrix.txt being a gene expression matrix with sample names and gene names.)\n";
}

open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
my $l = <IN>;

my @genes_q75 = ();
my @genes_q90 = ();
my @genes_q95 = ();

my $type = 1;
if ($ARGV[1] ne "") {
  $type = $ARGV[1];
}

while (my $l = <IN>) {

  my @a = split /\t/, $l, -1;

  my $n = shift @a;

  my $m = &median(\@a);
  my $s = &median_absolute_deviation(\@a);
  
  if ($s == 0) {
    next;
  }
  
  foreach my $r (@a) {
    $r = ($r - $m) / $s;
  }

  my @b = ();
  if ($type == 1) {
    @b = sort { $a <=> $b } @a;
  } else {
    @b = sort { $b <=> $a } @a;
  }

  my $q75 = int( 0.5 + 0.75 * @b );
  my $q90 = int( 0.5 + 0.90 * @b );
  my $q95 = int( 0.5 + 0.95 * @b );

  my @c_q75 = ($n, $b[$q75]);
  my @c_q90 = ($n, $b[$q90]);
  my @c_q95 = ($n, $b[$q95]);

  push @genes_q75, \@c_q75;
  push @genes_q90, \@c_q90;
  push @genes_q95, \@c_q95;

  #print "$n\t$b[$q75]\t$b[$q90]\t$b[$q95]\n";

}

if ($type == 1) {

  @genes_q75 = sort { $b->[1] <=> $a->[1] } @genes_q75;
  @genes_q90 = sort { $b->[1] <=> $a->[1] } @genes_q90;
  @genes_q95 = sort { $b->[1] <=> $a->[1] } @genes_q95;

} else {

  @genes_q75 = sort { $a->[1] <=> $b->[1] } @genes_q75;
  @genes_q90 = sort { $a->[1] <=> $b->[1] } @genes_q90;
  @genes_q95 = sort { $a->[1] <=> $b->[1] } @genes_q95;

}

my %COUNTS = ();
for (my $i=0; $i<100; $i++) {
  print STDERR "$genes_q95[$i]->[1]\n";
  $COUNTS{ $genes_q75[$i]->[0] } ++;
  $COUNTS{ $genes_q90[$i]->[0] } ++;
  $COUNTS{ $genes_q95[$i]->[0] } ++; 
}

foreach my $r (keys(%COUNTS)) {
  print "$r\t$COUNTS{$r}\n";
}



#
# MAD 
#
sub median_absolute_deviation {
  my ($a_ref) = @_;

  my $m = median($a_ref);    

  my @a = ();
  for (my $i=0; $i<@$a_ref; $i++) {
    my $t = abs( $a_ref->[$i] - $m );
    push @a, $t;
  }	  
  
  return median(\@a);
  
}

#
# compute the median value of an array
#
sub median {
    
  my ($a_ref) = @_;
  
  my $n     = scalar(@$a_ref);
  my @a_tmp = sort {$a <=> $b} @$a_ref;
  
  my $nm    = int($n / 2) - 1;  # floor - 1
    
  if ($n % 2 == 1)  {
    return $a_tmp[$nm + 1];
  } else {
    return ($a_tmp[$nm] + $a_tmp[$nm + 1]) / 2;
  }    

}
