use lib "$ENV{FIREDIR}/SCRIPTS";
use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

use Getopt::Long;
use Fire;
use Table;
use Sets;

use PostScript::Simple;
use strict;

my $scriptdir = "$ENV{FIREDIR}/SCRIPTS";

my $expfile     = undef;
my $data        = undef;
my $summaryfile = undef;
my $ps2pdf      = 1;
my $dologos     = 1;
my $h           = 30;
my $pvalues     = undef;
my $xsize       = undef;
my $ysep        = 100;
my $dna_rna     = 'DNA_RNA';

GetOptions ('expfile=s'              => \$expfile,
	    'data=s'                 => \$data,
	    'dologos=s'              => \$dologos,
	    'pvalues=s'              => \$pvalues,
	    'xsize=s'                => \$xsize,
	    'ysep=s'                 => \$ysep,
	    'dna_rna=s'              => \$dna_rna,
	    'summaryfile=s'          => \$summaryfile);



my $file = Sets::filename($expfile);

if (!defined($pvalues)) {
  $pvalues = "$expfile\_FIRE/$dna_rna/$file.matrix";
}


if (!defined($summaryfile)) {
  $summaryfile = "$expfile\_FIRE/$dna_rna/$file.summary";
}

my $outeps = "$expfile\_FIRE/$dna_rna/clusters.eps";
my $outpdf = "$expfile\_FIRE/$dna_rna/clusters.pdf";
my $outpng = "$expfile\_FIRE/$dna_rna/clusters.png";


my $ta = Table->new;

system("perl -pi -e 's/\\r//g' $expfile");


$ta->loadFile($expfile);
my $a_ref_clu = $ta->getArray();

my %CLUSTERS = ();
shift @$a_ref_clu;
foreach my $r (@$a_ref_clu) {
  push @{ $CLUSTERS{ $r->[1] } }, $r->[0];
}


my %PVALUES = undef;
my $nbclusters = undef;
if (defined($pvalues)) {
  print "loading $pvalues\n" if (-e $pvalues);
  my $h_ref_pv = $ta->getBidimensionalHash($pvalues);
  %PVALUES = %$h_ref_pv;

  $ta->loadFile($pvalues);
  $nbclusters = $ta->getNbColumns();
}


my %DATA = ();
$ta->loadFile($data);
my $a_ref_data = $ta->getArray();

my $A_REF_CONDS = shift @$a_ref_data;
shift @$A_REF_CONDS;
foreach my $r (@$a_ref_data) {
  my $n = shift @$r;
  $DATA{ $n } = $r;
}

#
# calculate the centroids
#
my %CENTROIDS = ();
foreach my $c (keys(%CLUSTERS)) {

  my $s        = [];
  my @cnt_rows = ();

  foreach my $r (@{$CLUSTERS{$c}}) {
    for (my $i=0; $i<@{ $DATA{$r} }; $i++) {
      if ($DATA{$r}->[$i] ne "") {
	$s->[$i] += $DATA{$r}->[$i];
	$cnt_rows[$i]++;
      } 
    }
    # $s = Sets::addArrays($s, $DATA{$r});    
  }

  my $n = @{$CLUSTERS{$c}};
  for (my $i=0; $i<@$s; $i++) { 
    
    if ($cnt_rows[$i] == 0) {
      $s->[$i] = 0.0;
    } else {
      $s->[$i] /= $cnt_rows[$i];
    }
  }
  $CENTROIDS{$c} = $s;

}

my $H_REF_SUM = Fire::loadFireMotifSummary($summaryfile);
my %CLUMOT    = ();
my $d = "TMP"; mkdir $d if (! -e $d);
my %LOGOS     = ();
my %NB2MOT    = ();
foreach my $m (keys(%$H_REF_SUM)) {
  my $clu = $H_REF_SUM->{$m}->{CLU};
  foreach my $c (@$clu) {
    push @{$CLUMOT{$c}}, $H_REF_SUM->{$m}->{CNT};
    print "$m is over-rep in cluster $c\n";
  }

  my $cnt = $H_REF_SUM->{$m}->{CNT};
  $NB2MOT{$cnt} = $m;

  if ($dologos == 1) {
    #
    # make motif logo
    #
    my $mo = Sets::myre2wm($m);
    if ($H_REF_SUM->{$m}->{RNA} == 1) {
      $mo =~ s/T/U/g;
    }
    
    
    open OUT, ">$d/$cnt.txt" or die "cannot open $d/$cnt.txt\n";
    print OUT $mo;
    close OUT;
    print "Outputing logo for $m\n";
    system("$scriptdir/weblogo/seqlogo -f $d/$cnt.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $d/$cnt.eps");
  }

  
  #
  # START: integrate Motif LOGO
  #  
  my $e  = new PostScript::Simple::EPS(file => "$d/$cnt.eps");

  # get height
  my $eh = $e->height;
  my $ew = $e->width;
 
  # height must be $h, so scale down to $h = k * LO
  $e->scale($h / $eh);

  $LOGOS{$cnt} = $e; 

}  


#
#  START DRAWING
#
#

my $xbase = 55;
my $ybase = 105;

if (!defined($xsize)) {
  $xsize = 800;
}


my $yclustsize = 150;
my $w          = 10;

if (!defined($nbclusters)) {
  $nbclusters = keys(%CENTROIDS);
}

my $ysize      = $nbclusters * ($yclustsize + $ysep); 

my $p          = new PostScript::Simple(xsize     => $xsize,
					ysize     => $ysize,
					colour    => 1,
					eps       => 1,
					units     => "pt");

$p->setlinewidth(0.5);


#
# draw cluster
#


$p->setcolour("black");
$p->setfont("Courier", 8);
my $i = 0;

my $myn = @{$CENTROIDS{0}};

$p->setfont("Arial", 14);
$p->text({ align => "center"}, $xbase + $myn * $w / 2,  $ysize - ($ybase - 30), "Cluster centroid");
$p->text({ align => "center"}, $xbase + $myn * $w + 100,  $ysize - ($ybase - 30), "Over-represented motifs");
$p->setfont("Courier", 8);

$ybase += 100;

foreach my $c (sort {$a <=> $b} (keys(%CENTROIDS))) {
  

  print "cluster '$c'\n";

  my $r = $CENTROIDS{$c};
  
  if (defined($CLUMOT{$c})) {
    
    my $ystart = $ybase + $i * ( $yclustsize + $ysep );


    $xbase -= 7;
    $p->setfont("Courier", 14);
    $p->text({ align => "right" }, $xbase - 10,  $ysize - ($ystart + 10), "C$c");
    $p->setfont("Courier", 8);
    
    $p->line($xbase, $ysize - ($ystart + 10), $xbase, $ysize - ($ystart + $yclustsize - 10));
 
   $p->line($xbase - 3, $ysize - ($ystart +     $yclustsize / 4 ), $xbase + 3, $ysize - ($ystart +     $yclustsize / 4));
   $p->line($xbase - 3, $ysize - ($ystart + 3 * $yclustsize / 4 ), $xbase + 3, $ysize - ($ystart + 3 * $yclustsize / 4));
   $p->line($xbase - 3, $ysize - ($ystart + 2 * $yclustsize / 4 ), $xbase + 3, $ysize - ($ystart + 2 * $yclustsize / 4));
    
    $p->text({ align => "right" }, $xbase - 3, $ysize - ($ystart +     $yclustsize / 4 + 3),  "1.0");
    $p->text({ align => "right" }, $xbase - 3, $ysize - ($ystart + 3 * $yclustsize / 4 + 3), "-1.0"); 
    $p->text({ align => "right" }, $xbase - 3, $ysize - ($ystart + 2 * $yclustsize / 4 + 3),  "0.0");
    
    

    $xbase += 7;
    
    
    my $n = @{$CENTROIDS{$c}};
    for (my $j=0; $j<$n; $j++) {
      my $xi = $xbase + $j * $w;
      my $yi = $ysize - ( $ystart + $yclustsize / 2);
      my $xj = $xi+$w;
      my $yj = $yi + $r->[$j] * $yclustsize / (2 * 2);
      #print "$r->[$j]\n";

      #print "$xi, $yi, $xj, $yj\n";
      
      $p->box($xi, $yi, $xj, $yj);
      
      my @a = split /\_/, $A_REF_CONDS->[$j];
      #shift @a; shift @a; 
      #pop @a; 
      my $t = $a[0];

      
      $t =~ s/Biological\ Replicate//g;
      $t =~ s/Repeated\ Hybridization/RepHyb/g;
      $t =~ s/\ +/\ /g;

      $p->text({rotate => 90}, $xi + 2*$w/3, Sets::max($yi, $yj) + 2, $t);
      
    }
    #<STDIN>;
    my $j = 0;
    foreach my $m (@{$CLUMOT{$c}}) {
      $p->_add_eps($LOGOS{$m}, $xbase + $n*$w + 10 + $j * 50,  $ysize - ( $ystart + $yclustsize / 2 + 10)); 
      if (defined($pvalues)) {
	my $mm = $NB2MOT{$m}; $mm =~ s/U/T/g;
	print "PVALUES{$mm}{C$c} = " . $PVALUES{$mm}{"C$c"} . "\n";
	$p->text($xbase + $n*$w + 20 + $j * 50,  $ysize - ( $ystart + $yclustsize / 2 + 10), "p<1e-" . int($PVALUES{$mm}{"C$c"}));     
      }
      $j++;
    }

    
    
    $i++;
  }
  
  
  #last if ($i == 5);
}

$p->output($outeps);


if ($ps2pdf == 1) {
  print "Creating PDF $outpdf  ... ";
  system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf");
  system("pstoimg -antialias $outeps");
  print "Done.\n";
}


