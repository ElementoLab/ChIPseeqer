use lib "$ENV{FIREDIR}/SCRIPTS";

use Sets;
use Table;
use strict;
use Fasta;
use strict;

if (@ARGV == 0) {
  die "Usage: perl quantize_expression_vector.pl -expfile -outfile -mbins -ebins\n";
}

my $expfile         = Sets::get_parameter(\@ARGV, "-expfile");


my $divbins      = 50.0;
if (Sets::exist_parameter(\@ARGV, "-divbins") == 1) {
  $divbins       = Sets::get_parameter(\@ARGV, "-divbins");
}

my $outfile         = undef;
if (Sets::exist_parameter(\@ARGV, "-outfile") == 1) {
  $outfile         = Sets::get_parameter(\@ARGV, "-outfile");
}

my $mbins      = 2;
if (Sets::exist_parameter(\@ARGV, "-mbins") == 1) {
  $mbins       = Sets::get_parameter(\@ARGV, "-mbins");
}

my $ebins      = undef;
if (Sets::exist_parameter(\@ARGV, "-ebins") == 1) {
  $ebins       = Sets::get_parameter(\@ARGV, "-ebins");
}



open OUT, ">$outfile";


my %expgenes_good = ();
open IN, $expfile;
my $l = <IN>;
print OUT $l;

my @INGENES = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;
  $expgenes_good{ $a[0] } = $a[1];
  push @INGENES, $a[0];
}
close IN;

#
# do a loop while there are differences between iterations
#

my $a_ref_fil               = [];
my $a_ref_fil_new           = [];
my $h_ref_expgenes_good     = \%expgenes_good;
my $h_ref_expgenes_good_new = undef;
  
my $N = scalar( keys(%$h_ref_expgenes_good) );

# how many bins ?

if (!defined($ebins)) {
  $ebins = int( 0.5 + $N / ( $divbins * $mbins ) ); 
}

my @VAL = ();
my @GEN = ();
foreach my $g (sort(keys(%$h_ref_expgenes_good))) {
  push @VAL, $h_ref_expgenes_good->{$g};
  push @GEN, $g;
}

my $vq = Quantize(\@VAL, $ebins);

my %H  = ();
for (my $i=0; $i<$N; $i++) {
  $H{ $GEN[$i] } = $vq->[$i];
}


foreach my $g (@INGENES) {
  if (defined($H{$g})) {
    print OUT "$g\t$H{$g}\n";
  }
}

close OUT;


sub Quantize{ 

  my ($v, $D) = @_;
  
  
  my $n = scalar( @$v );
  
  my $binsize = int( 0.5 + $n / $D);

  my @fi = ();
  for (my $i=0; $i<$n; $i++) {
    $fi[$i]->[0] = $i;
    $fi[$i]->[1] = $v->[$i];
  }
  
  @fi = sort { $a->[1] <=> $b->[1] } @fi;

  my $qv = [];

  for (my $i=0; $i<$D; $i++) {
    my $j        = $i * $binsize;

    while (( $j < ($i+1)*$binsize) && ($j < $n)) {
      $qv->[ $fi[$j]->[0] ] = $i;
      $j ++;
    }
  }

  return $qv;
}


