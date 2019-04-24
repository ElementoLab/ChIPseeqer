#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Matrix2EPS;
use AggloClust;

use strict;

#
# read in the interaction matrix
#

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my @MATRIX = ();
my %CONDS  = ();
my $numconds = 0;

my @names  = ();

foreach my $r (@$a_ref) {

  #print "$r->[0] / $r->[1]\n";

  my $idx1 = undef;
  if (defined($CONDS{$r->[0]})) {
    $idx1 = $CONDS{$r->[0]};
  } else {
    $idx1 = $numconds++;
    $CONDS{$r->[0]} = $idx1;
    $names[$idx1] = $r->[0];
  }

  my $idx2 = undef;
  if (defined($CONDS{$r->[1]})) {
    $idx2 = $CONDS{$r->[1]};
  } else {
    $idx2 = $numconds++;
    $CONDS{$r->[1]} = $idx2;
    $names[$idx2] = $r->[1];
  }

  my $z    = $r->[4] * Sets::sign($r->[5]);

  $MATRIX[ $idx1 ][ $idx2 ] = $z;
  $MATRIX[ $idx2 ][ $idx1 ] = $z;



}

for (my $i=0; $i<$numconds; $i++) {
  $MATRIX[ $i ][ $i ] = 1000;
}

#
# cluster profiles based on z-scores
#
my $ac = AggloClust->new;
$ac->setAlgoClust("min");

$ac->setUseAbs(1);
$ac->setUseCorr(1);

$ac->setDistanceMatrix(\@MATRIX);
$ac->agglomerate_using_max_linkage();

my $a_ref_o = $ac->getDFSOrder();

print join("-", @$a_ref_o) . "\n";

my $n          = @MATRIX;
my $a_ref_newm = [];
my $a_ref_newn = [];
my $a_ref_newo = [];


for (my $j=0; $j<@$a_ref_o; $j++) {
  
  # reorganize row  $a_ref_o->[ $j ]
  my $r = [];
  for (my $k=0; $k<@$a_ref_o; $k++) {
    $r->[ $k ] = $MATRIX[ $a_ref_o->[ $j ] ][ $a_ref_o->[ $k ] ];
  }

  # send it where it belongs
  $a_ref_newm->[ $j ] = $r;

  $a_ref_newn->[ $j ] = $names[ $a_ref_o->[ $j ] ];
  $a_ref_newo->[ $j ] =  $a_ref_newn->[ $j ];
}


#
# show it
#
my $ma = Matrix2EPS->new;
$ma->setOutputFileName("test.eps");
$ma->setMin(-50);
$ma->setMax(50);
$ma->setMinText("Down");
$ma->setMaxText("Up");
$ma->setVerbose(1);
$ma->setH(10);
$ma->setW(10);

$ma->setColMap("/Users/olivier/PROGRAMS/MYTREEVIEW/HEATMAPS/rbg_cmap.csv");
$ma->drawRowNames("right");

$ma->setMatrix($a_ref_newm, $a_ref_newn, $a_ref_newo);

$ma->setXbase(50);
$ma->setYbase(300);
$ma->setXright(300);

$ma->setRownameFontSize(10);
$ma->setHeaderFontSize(10);


$ma->draw();

$ma->drawScale(20,280);

$ma->output();
$ma->pdfify();

