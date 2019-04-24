use lib "$ENV{FIREDIR}/SCRIPTS";
use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

use Table;
use Sets;
use Getopt::Long;
use PostScript::Simple;
use strict;

my $progdir        = "$ENV{FIREDIR}/PROGRAMS";

if (@ARGV == 0) {
  die "Usage : perl mi_create_dist_pvalues.pl -expfile FILE -distrepfile FILE -summaryfile FILE -profiles FILE -motifdir DIR -quantized INT -outdistmatrix FILE -outoriematrix FILE\n";
}

my $distrepfile = Sets::get_parameter(\@ARGV, "-distrepfile");
my $expfile     = Sets::get_parameter(\@ARGV, "-expfile");

my $mapdir = "$distrepfile\_OUT";
if (! -e $mapdir) {
  mkdir $mapdir;
}


my $profiles      = Sets::get_parameter(\@ARGV, "-profiles");

my $summaryfile   = Sets::get_parameter(\@ARGV, "-summaryfile");

my $outdistmatrix = Sets::get_parameter(\@ARGV, "-outdistmatrix");
my $outoriematrix = Sets::get_parameter(\@ARGV, "-outoriematrix");


my $motifdir    = undef;
if (Sets::exist_parameter(\@ARGV, "-motifdir") == 1) {
  $motifdir    = Sets::get_parameter(\@ARGV, "-motifdir");
}

my $quantized   = Sets::get_parameter(\@ARGV, "-quantized");
my $ps2pdf      = 1;
my $rootdir     = ".";
my $seqlen      = Sets::get_parameter(\@ARGV, "-seqlen");


my $rightlabel  = 'right';
if (Sets::exist_parameter(\@ARGV, "-rightlabel") == 1) {
  $rightlabel    = Sets::get_parameter(\@ARGV, "-rightlabel");
}

my $leftlabel   = 'left';
if (Sets::exist_parameter(\@ARGV, "-leftlabel") == 1) {
  $leftlabel    = Sets::get_parameter(\@ARGV, "-leftlabel");
}

my $limit       = undef;




# traverse all motifs
#   if they have a distance bias, output matrix
#   if pos / ori bias, output motif map


my $ta = Table->new;


#
#  read in expression data
#
$ta->loadFile($expfile);
my $h_ref_exp = $ta->getIndexKV(0,1);

#
#  end read in exp data
#

#
#  read in the profiles
#
$ta->loadFile($profiles);
my $a_ref_prof = $ta->getArray();
#
#  end read in profiles
#


#
#  read in the summary file
#

$ta->loadFile($summaryfile);
my $a_ref_mo = $ta->getArray();
my %STAT      = ();
my @MOTIFS    = ();

foreach my $r (@$a_ref_mo) {

  push @MOTIFS, $r->[0];
  
  my %a_tmp = ( "RNA"    => $r->[1],
		"COPIES" => $r->[2],
		"MI"     => $r->[3],
		"RANK"   => $r->[4], 
		"Z"      => $r->[5], 
		"R"      => $r->[6], 
		"S"      => $r->[7],
		"SEED"   => $r->[8],
		"DIST"   => $r->[9],
		"ORIE"   => $r->[10],
		"CONS"   => $r->[11] );
  
  #
  # enriched clusters
  #
  my @a_clu = ();
  for (my $i=12; $i<@$r; $i++) {
    push @a_clu, $r->[$i];
  }

  $a_tmp { "CLUSTERS" } = \@a_clu;

  $STAT{ $r->[0] } = \%a_tmp;   
  
}


#  
# load distance report file  
#


print "Parsing distance report file ... ";

$ta->loadFile($distrepfile);
my $a_ref = $ta->getArray();


#
#  load the distance data (mbins)
#
my %DISTDATA = ();
my %MBINS    = ();
my %ORIEDATA = ();
my %EBINS    = ();
foreach my $r (@$a_ref) {
  
  my $re = shift @$r;
  my $ty = shift @$r;
  my $bi = shift @$r;

  
  #
  # distance 
  #
  if ($ty eq "d_avg") {
    if ($bi eq "nan") {
      $MBINS{$re} = $r;
    } else {
      my $mm = @$r;
      if (($quantized == 0) || ($r->[$mm-1] =~ /[\]\[]/)) {
	$EBINS{$bi} = pop @$r;
      }
      $DISTDATA{ $re }[$bi] = $r;
    }
  }

  if ($ty eq "o5") {
    $ORIEDATA{$re}[$bi][0] = $r->[1]; 
  }
  
  if ($ty eq "o3") {
    $ORIEDATA{$re}[$bi][1] = $r->[1]; 
  }
}

print " Done.\n";


#
# outdistmatrix, outoriematrix
#
open OD, ">$outdistmatrix" or die "Cannot open $outdistmatrix.\n";
open OO, ">$outoriematrix" or die "Cannot open $outoriematrix.\n";


#
# traverse all motifs
#

my $cnt = 0;
foreach my $re (@MOTIFS) {


  #
  #  IMPORTANT : SKIP MOTIFS WITH WITH NO ORIENTATION OR DISTANCE BIASES
  #
  if (($STAT{ $re }->{DIST} != 1) && ($STAT{ $re }->{ORIE} == 0)) {
    print "$re has no position or orientation bias.\n";
    $cnt ++; next;
  }

  
  if ($STAT{$re}->{ORIE} >= 1) {

    my $s = $ORIEDATA{$re};
    my $n = @$s; # number of clusters
    

    # calculate the basal frequencies
    my $f1 = 0; my $f2 = 0;

    foreach my $t (@$s) {
      $f1 += $t->[0];
      $f2 += $t->[1];
    }

    my $sum = $f1 + $f2;
    $f1 = $f1 / $sum;  $f1 = 0.5;
    $f2 = $f2 / $sum;  $f2 = 0.5;

    my $cnt = 0;
    foreach my $t (@$s) {

      my $sum = $t->[0]+$t->[1];
      my $p1m = undef;
      my $p1  = binom_test_greater( $t->[0], $sum, $f1, \$p1m);
      my $p2m = undef;
      my $p2  = binom_test_greater( $t->[1], $sum, $f2, \$p2m);

      if ($p1 < $p1m) {  # if (p1 is greater then p1m, we have an enrichment)
	$p1 = sprintf("%4.3f", - Sets::log10($p1));
      } else {
	$p1 = sprintf("%4.3f", + Sets::log10($p1m));
      }

      if ($p2 < $p2m) {  # if (p1 is greater then p1m, we have an enrichment)
	$p2 = sprintf("%4.3f", - Sets::log10($p2));
      } else {
	$p2 = sprintf("%4.3f", + Sets::log10($p2m));
      }

      print OO "$re\t$cnt\t$p1\t$p2\n";

      $cnt ++;

    }
  }  # end if orientation bias


  #
  # DISTANCE count matrix
  #
  my $r = $DISTDATA{$re};
  my $n = @$r; 
  my $m = @{ $r->[0] };

  if ($STAT{$re}->{DIST} == 1) {

    print "Processing matrix $n x $m.\n";

    my @SUMS_N = ();
    my @SUMS_M = ();
    
    my $j = 0;
    foreach my $s (@$r) {
      $SUMS_N[$j] = Sets::arraySum($s);
      $j ++;    
      for(my $i=0; $i<@$s; $i++) {
	$SUMS_M[$i] += $s->[$i];
      }    
    }
    
    my $N = Sets::arraySum(\@SUMS_N);
    
    for(my $i=0; $i<@SUMS_N; $i++) {
      $SUMS_N[$i] /= $N;
    }
    
    for(my $i=0; $i<@SUMS_M; $i++) {
      $SUMS_M[$i] /= $N;
    }
  
    
    my @M = ();
    for(my $i=0; $i<@SUMS_N; $i++) {
      for(my $j=0; $j<@SUMS_M; $j++) {

	
	my $o  = $N * $SUMS_N[$i] * $SUMS_M[$j];
	my $bk = $r->[$i]->[$j];
	my $bn = $N;
	my $bp = $SUMS_N[$i] * $SUMS_M[$j];
	my $p2 = undef;
	my $p1 = binom_test_greater( $bk, $bn, $bp, \$p2);
	
	my $p  = undef;
	
	if ($p1 < 1e-100) {
	  $p1 = 1e-100;
	}
	if ($p2 < 1e-100) {
	  $p2 = 1e-100;
	}
	if ($p1 < $p2) {
	  $p = sprintf("%4.3f", - Sets::log10($p1));
	} else {
	  $p = sprintf("%4.3f", + Sets::log10($p2));
	}
	
	$M[$i][$j] = $p;

	#if ($p > - Sets::log10( 0.05 / $n)) {
	
	print OD "$re\t$i\t" . $MBINS{$re}->[$j] . "\t$p\n"; 
	#}
	
      }
    }

  } # end if dist bias   

  $cnt ++;
  
  last if (defined($limit) && ($cnt == $limit));
}



close OD;
close OO;



sub binom_test_greater {
  
  my ($k, $n, $p, $p1) = @_;

  if (($n == 0) && ($k == 0)) {
    $$p1 = 1.0;
    return 1.0;
  }

  my $todo = "$progdir/binom_test_greater $k $n $p 1";
  my $out = `$todo`;  
  $out =~ s/[\n\r]//g;

  my @a = split /\t/, $out, -1;
  
  $out = $a[0];
  $$p1 = $a[1];
  
  return $out;
}
