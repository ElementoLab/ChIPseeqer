# takes as input summary file (with over-rep clusters)
#                signif rep file (with fractions)

use Fire;
use Sets;
use strict;
use Data::dumper;

my $h_ref_s = Fire::loadFireMotifSummary($ARGV[0]);
my $h_ref_r = Fire::loadFireMotifRep($ARGV[1]);

foreach my $m (keys(%$h_ref_s)) {

  #print "M=$m\n";
  
  my @clu = @{$h_ref_s->{$m}->{CLU}};
  my $k = @clu;
  # print "present in " . scalar(@clu) . "\n";
  
  # get number of clusters

 

  #print Dumper($h_ref_r->{$m});

  my $n = @{ $h_ref_r->{$m}->{F} };

  my $tot = 0;
  my $in  = 0;
  my $totin = 0;
  my $out = 0;
  my $totout = 0;

  for (my $i=0; $i<$n; $i++) {
    my $p = $h_ref_r->{$m}->{F}->[$i] * $h_ref_r->{$m}->{N}->[$i];
    if (Sets::in_array($i, @clu)) {
      $in += $p;
      $totin += $h_ref_r->{$m}->{N}->[$i];
    } else {
      $out += $p;
      $totout += $h_ref_r->{$m}->{N}->[$i];

    }
    $tot += $p;
  }

  if ($totin > 0) {
    $in  = sprintf("%3.2f", $in /$totin);
    $out = sprintf("%3.2f", $out/$totout);
    print "$m\t$k\t$in\t$out\n";
  }

}

