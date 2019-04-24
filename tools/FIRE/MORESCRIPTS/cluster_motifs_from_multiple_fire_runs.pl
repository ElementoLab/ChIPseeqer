if (!$ENV{FIREDIR} || ($ENV{FIREDIR} eq '')) {
  die "Please define $ENV{FIREDIR}.\n";
}

use lib "$ENV{FIREDIR}/SCRIPTS";

use Sets;
use Table;
use Getopt::Long;
use AggloClust;
use strict;

my $expfiles = undef;
my $rna      = undef;
my $tc       = 0.7;
my $showall  = 1;

my $outfile  = undef;
my $outfull  = undef;
my $listexpfiles = undef;
my $verbose  = 1;
my $fastafile = undef;
my $doprofiles=1;

if (@ARGV == 0) {
  die "Usage: perl .. --expfiles --tc --rna=0/1 --showall\n";
}

GetOptions ('expfiles=s'             => \$expfiles,
	    'listexpfiles=s'         => \$listexpfiles,
	    'outfile=s'              => \$outfile,
	    'outfull=s'              => \$outfull,
	    'tc=s'                   => \$tc,
	    'fastafile=s'            => \$fastafile,
	    'doprofiles=s'           => \$doprofiles,
	    'rna=s'                  => \$rna,
	    'showall=s'              => \$showall);


my $ta     = Table->new;

#1 build list from .summary files
my $a_ref_f = undef;
if (defined($expfiles)) {
  $a_ref_f = Sets::getFiles($expfiles);
} elsif (defined($listexpfiles)) {
  $ta->loadFile($listexpfiles);
  $a_ref_f = $ta->getColumn(0);
}




my %SEEDS  = ();

foreach my $f (@$a_ref_f) {
  my $fi = Sets::filename($f);
  my $ff = "$f\_FIRE/DNA_RNA/$fi.summary";
  
  if (-e $ff) {
    $ta->loadFile($ff);

    if ($verbose == 1) {
      print "Loading $ff.\n";
    }
    
    my $a_ref = $ta->getArray();
    
    foreach my $r (@$a_ref) {
      if (defined($r->[12]) && ($r->[12] ne "")) {	
	
	next if ($r->[1] != $rna);
	
	# only add to the list if seed was not seen before, or if z > z of previous seed
	if (!defined($SEEDS{$r->[8]}) || (defined($SEEDS{$r->[8]}) && ($r->[5] > $SEEDS{$r->[8]}->{Z}))) { 
	  my %H = ("RE" => $r->[0], "Z" => $r->[5]);
	  $SEEDS{ $r->[8] } = \%H;
	}
	
      }
    }
    
  }
}

my @MOTIFS = values(%SEEDS);

my $txt = "";
my $cnt = 1;
foreach my $m (@MOTIFS) {
  my $ar = Sets::get_array_from_re($m->{RE});
  my $nr = @$ar;
  my $wm = Sets::myre2wm($m->{RE});
  $txt .= "Motif $cnt\n$wm"; $txt .= "*" x $nr; $txt .= "\n\n";
  $cnt ++;
}


print "Found " . scalar(@MOTIFS) . " motifs\n";



my $tmpfile1 = Sets::getTempFile("tmp-cluster.ace");
my $tmpfile2 = Sets::getTempFile("tmp-cluster-compareace.txt");

open OUT, ">$tmpfile1";
print OUT $txt;
close OUT;

my $todo = "$ENV{FIREDIR}/PROGRAMS/MyCompareAce $tmpfile1 $tmpfile1";
if ($rna == 1) {
  $todo .= " -ss ";
}
$todo .= " > $tmpfile2";
system($todo);



# reading matrix in
$ta->loadFile($tmpfile2);
my @DIST = ();
my $a_ref = $ta->getArray();
foreach my $r (@$a_ref) {

  # $DIST[ $r->[1] - 1 ][ $r->[3] - 1 ] = 1 - $r->[4];
  # $DIST[ $r->[3] - 1 ][ $r->[1] - 1 ] = 1 - $r->[4];

  $DIST[ $r->[0] ][ $r->[1] ] = 1 - $r->[2];
  $DIST[ $r->[1] ][ $r->[0] ] = 1 - $r->[2];

  
}

my $ac = AggloClust->new;
$ac->setDistanceMatrix(\@DIST);
$ac->setMinD(1-$tc);

my $a_ref_c = $ac->agglomerate_using_avg_linkage();


open OUT1, ">$outfile" or die "Cannot open $outfile.\n";
if (defined($outfull)) {
  open OUT2, ">$outfull" or die "Cannot open $outfull.\n";
}

my @MOTIFS_REP   = ();
my @MOTIFS_REP_M = ();
my %MOTIFS_REP_Z = ();
my %MOTIFS_REP_KEEP = ();

foreach my $r (@$a_ref_c) {

  my $b_z = -1; 
  my $b   = undef;
  foreach my $s (@$r) {
    if ($MOTIFS[$s]->{"Z"} > $b_z) {
      $b_z = $MOTIFS[$s]->{"Z"};
      $b   = $s;
    }	
  }

  push @MOTIFS_REP, $b;
  push @MOTIFS_REP_M, $MOTIFS[$b]->{RE};
  $MOTIFS_REP_Z{ $MOTIFS[$b]->{RE} } = $MOTIFS[$b]->{Z};
  $MOTIFS_REP_KEEP{ $MOTIFS[$b]->{RE} } = 1;

  print OUT1 $MOTIFS[$b]->{"RE"};
  print OUT2 $MOTIFS[$b]->{"RE"};

  if ($showall == 1) {
    foreach my $s (@$r) {
      print OUT2 "\t" . $MOTIFS[$s]->{"RE"} if ($s != $b);
    }
  }

  print OUT1 "\n";
  print OUT2 "\n";
}

close OUT1;
close OUT2;

unlink $tmpfile1;
unlink $tmpfile2;


# clustering, level two
#  look at all pairs in @MOTIFS_REP .. in coccur, remove weakest (smallest z-score)



# if fastafile defined, create .profiles file ...
if (defined($fastafile)) {

  print "Generating profiles .. \n";
  my $tmpfile0 = "tmp-combined-motifs.txt";
  Sets::writeSet(\@MOTIFS_REP_M, $tmpfile0);

  my $todo = "perl $ENV{FIREDIR}/SCRIPTS/generate_motif_profiles.pl -fastafile $fastafile -summaryfile $tmpfile0 -rna $rna -noflank 1 -outfile $tmpfile0.profiles ";
  
  if ($doprofiles == 1) {
    system("$todo");
  }

  $ta->loadFile("$tmpfile0.profiles");
  my $a_ref = $ta->getArray();

  my %POS = ();
  foreach my $r (@$a_ref) {
    push @{ $POS{ $r->[0] }->{ $r->[1] } }, $r->[2];
  }

  my @motifs = keys(%POS);

  for (my $i=0; $i<@motifs-1; $i++) {
    
    my $tmp1       = Sets::get_array_from_re($motifs[$i]);
    my $k1         = @$tmp1;
    my $h_ref_gp1 = $POS{$motifs[$i]};
    my @g1        = keys( %$h_ref_gp1 );
    my $n1        = 0;

    foreach my $h (@g1) {
      $n1 += @{$h_ref_gp1->{$h}};
    }
    
    
    #print "num genes in g1 $h_ref_gp1 = " . scalar(@g1) . "\n";
    
    for (my $j=$i+1; $j<@motifs; $j++) {
      
      my $n2        = 0;
      
      my $tmp2       = Sets::get_array_from_re($motifs[$j]);
      my $k2         = @$tmp2;
      
      my $h_ref_gp2 = $POS{$motifs[$j]};
      my @g2        = keys( %$h_ref_gp2 );

      foreach my $h (@g2) {
	$n2 += @{$h_ref_gp2->{$h}};
      }

      my $gg        = Sets::getOverlapSet(\@g1, \@g2);

      #print "Overlap has " . scalar(@$gg) . "\n";

      my $ov        = 0;

      foreach my $g (@$gg) {
	foreach my $p1 (@{$h_ref_gp1->{$g}}) {
	  foreach my $p2 (@{$h_ref_gp2->{$g}}) {
	    if (Sets::sequencesOverlap($p1, $p1+$k1, $p2, $p2+$k2)) {
	      $ov ++;
	    }
	  }
	}
      }

      if (($ov/$n1>0.33) || ($ov/$n2>0.33)) {

	#print "$motifs[$i] and $motifs[$j]: ov=$ov, n1=$n1, n2=$n2\n";

	if ($MOTIFS_REP_Z{ $motifs[$i] } < $MOTIFS_REP_Z{ $motifs[$j] }) {
	  $MOTIFS_REP_KEEP{ $motifs[$i] } = 0;
	} else {
	  $MOTIFS_REP_KEEP{ $motifs[$j] } = 0;
	}

      }

    }

  }


  open OUTO, ">$outfile.nooverlap";
  foreach my $m (keys(%MOTIFS_REP_KEEP)) {
    print OUTO "$m\n" if ($MOTIFS_REP_KEEP{$m} == 1);
  }
  close OUTO;
  
}

