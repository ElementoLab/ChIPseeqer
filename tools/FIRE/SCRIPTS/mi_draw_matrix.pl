use lib "$ENV{FIREDIR}/SCRIPTS";
use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

use Table;
use Sets;
use Getopt::Long;
use PostScript::Simple;
use strict;


my $scriptdir      = "$ENV{FIREDIR}/SCRIPTS";
my $progdir        = "$ENV{FIREDIR}/PROGRAMS";

my $matfile        = undef;
my $densityfile    = undef;
my $motiffile      = undef;
my $quantized      = 0;
my $clustfile      = undef;
my $rna            = 0;
my $sortcolsbysize = 0;
my $sortrowsbymax  = 0;
my $columnsfile    = undef;
my $distfile       = undef;
my $limit          = undef;
my $ps2pdf         = 0;
my $w              = undef;
my $h              = undef;
my $every          = 1;
my $consfile       = undef;
my $gofile         = undef;
my $colmap         = undef;
my $dopng          = 0;
my $rootdir        = '.';;
my $showstatlabels = 1;
my $showscale      = 1;
my $ybase          = undef;
my $redoweblogo    = 1;
my $motifnames     = undef;
my $outeps         = undef;
my $lp_t_draw      = 20;
my $expfile        = undef;

if (@ARGV == 0) {
  die "Usage: perl mi_draw_matrix.pl  --matfile=FILE --symmaryfile=FILE --clustfile=FILE --columnsfile=FILE --gofile=FILE --ps2pdf=INT --every=1 --quantized=INT\n";
}

GetOptions (
	    'expfile=s'        => \$expfile,
	    'matfile=s'        => \$matfile,
	    'densityfile=s'    => \$densityfile,
	    'outeps=s'         => \$outeps,
	    'summaryfile=s'    => \$motiffile,
	    'clustfile=s'      => \$clustfile,
	    'motifnames=s'     => \$motifnames,
	    'columnsfile=s'    => \$columnsfile,
	    'distfile=s'       => \$distfile,
	    'consfile=s'       => \$consfile,
	    'sortcolsbysize=s' => \$sortcolsbysize,
	    'sortrowsbymax=s'  => \$sortrowsbymax,
	    'showstatlabels=s' => \$showstatlabels,
	    'showscale=s'      => \$showscale,
	    'w=s'              => \$w,
	    'h=s'              => \$h,
	    'gofile=s'         => \$gofile,
	    'every=s'          => \$every,
	    'limit=s'          => \$limit,
	    'ps2pdf=s'         => \$ps2pdf,
	    'colmap=s'         => \$colmap,
	    'dopng=s'          => \$dopng,
	    'ybase=s'          => \$ybase,
	    'rootdir=s'        => \$rootdir,
	    'redoweblogo=s'    => \$redoweblogo,
	    'quantized=s'      => \$quantized,
	    'lp_t_draw=s'      => \$lp_t_draw);


if (! -e $clustfile) {
  print "$clustfile does not exist, ignoring.\n";
  $clustfile = undef;
}

my $ta = Table->new;


#
#  load GO (for clusters)
#
my %GO = ();
if (defined($gofile)) {
  $ta->loadFile($gofile);
  %GO = %{ $ta->getIndexKV(0,1) };
  print "GO file loaded.\n";
}
  

my @MOTIFS = ();
my $a_ref_clust = undef;


#
#  read in the matrix file
#
$ta->loadFile($matfile);
my $a_ref_M      = $ta->getArray();
if (!defined($a_ref_clust)) {
  @MOTIFS = @{ $ta->getColumn(0) };
  shift @MOTIFS;
}
my $a_ref_header = shift @$a_ref_M; shift @$a_ref_header;
my %MATRIX       = ();
foreach my $r (@$a_ref_M) {
  my $m = shift @$r;
  $MATRIX{ $m } = $r;
}



#
#  read in the density file
#
my %DENSITIES       = ();
if (defined($densityfile)) {
  my $a_ref_D = undef;
  $ta->loadFile($densityfile);
  $a_ref_D      = $ta->getArray();
  shift @$a_ref_D; 
  foreach my $r (@$a_ref_D) {
    my $m = shift @$r;
    $DENSITIES{ $m } = $r;
  }
}


#
# read in motifs names
#
$ta->loadFile($motifnames);
my $h_ref_names   = $ta->getIndexKV(0,1);


#
# read in the clusters
#
if (defined($clustfile)) {
  $ta->loadFile($clustfile);
  $a_ref_clust = $ta->getColumn(1);
  my $col      = $ta->getColumn(0);
  @MOTIFS      = @$col;
}


#
#  read in the summary file
#
$ta->loadFile($motiffile);
my $a_ref_mo = $ta->getArray();
my %STAT         = ();
my %MOTIF_NUMBER = ();

my $cnt = 0;
foreach my $r (@$a_ref_mo) {
 
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
		"CONS"   => $r->[11],
		"NAME"   => $h_ref_names->{ $r->[0] });

  $STAT{ $r->[0] }         = \%a_tmp;   
  $MOTIF_NUMBER{ $r->[0] } = $cnt ++;
}



$ta->loadFile($columnsfile);
my $a_ref_cols  = $ta->getColumn(0);
for (my $i=0; $i<@$a_ref_cols; $i++) {
  $a_ref_header->[$i] .= " " . $GO{ $a_ref_cols->[$i] };
  if (defined($GO{ $a_ref_cols->[$i] })) {
    #print "GO annotation added for cluster $a_ref_cols->[$i].\n";
  } else {
    #print "No GO annotation for $a_ref_cols->[$i].\n";
  }	
}




#
# load color map
#
my $A_REF_COLMAP = undef;
if (defined($colmap)) {
  
  $ta->setDelim(" ");
  $ta->loadFile($colmap);
  $A_REF_COLMAP = $ta->getArray();
  $ta->setDelim('\t');
}


#$a_ref_M = sort2DArrayByMax($a_ref_M);



print "Now doing the graphical display.\n";

#
#  START DRAWING
#
#

die "No motifs in \@MOTIFS ...\n" if (@MOTIFS == 0);


my $xbase        = 35;
if (!defined($ybase)) {
  $ybase        = 250;
}

if (!defined($h)) {
  $h = 40;
}
my $ysize        = $ybase + $h * scalar(@MOTIFS) + 50;
$ysize = Sets::max($ysize, 370);
my $xsize        = 850;

print "xsize = $xsize, ysize = $ysize, xbase = $xbase, ybase = $ybase\n";

if (!defined($w)) {
  $w = int( 0.5 + ((2 * $xsize / 5.0) + 5)  / @$a_ref_header );
} 

my $p = new PostScript::Simple(xsize => $xsize, #papersize => "A4",
			       ysize => $ysize, #papersize => "A4",
			       #$papersize => "A4",
			       colour    => 1,
			       #landscape    => 1,
			       #coordorigin => "LeftTop",
			       #direction   => "RightDown",
			       eps       => 1,
			       units     => "pt");

$p->setlinewidth(0.5);



my $d = "$motiffile" . "_OUT";
if (! -e $d) {
  mkdir $d;
}



#
# show header
#
my $j = 0; 
my $n = @$a_ref_header;
my $gheader = 1;
if (($quantized == 1) || ($gheader == 0)) {
  foreach my $c (@$a_ref_header) {
    
    if ((($j % $every == 0) && ($j < $n-$every)) || ($j == $n-1)) {
      $p->setcolour("black");
      $p->setfont("Courrier", 6);
      $p->text( { rotate => 90 }, $xbase + $j * $w + 3*$w/4, $ysize - ($ybase-3), $c);  
      
    }
    $j ++;
  }
} else {

  # graphical header for continuous data
  
  my $min_i1 =  1000000;
  my $max_i2 = -1000000;
  my @bins   = ();
  foreach my $c (@$a_ref_header) {
    my ($i1, $i2) = $c =~ /\[(.+?)\;(.+?)\]/;    
    $min_i1 = $i1 if ($i1 < $min_i1);
    $max_i2 = $i2 if ($i2 > $max_i2);
    my @a_tmp = ($i1, $i2); push @bins, \@a_tmp;
  }
  
  my $th = $max_i2 - $min_i1;
  my $hi = $h * 1.5;

  $j = 0;
  foreach my $c (@$a_ref_header) {
    
    my $h1 = $hi * ($bins[$j]->[0] - $min_i1 ) / $th ;
    my $h2 = $hi * ($bins[$j]->[1] - $min_i1 ) / $th ;

    $p->setcolour("black");    
    $p->box({filled => 1}, 
	    $xbase + $j * $w,      $ysize - ($ybase - 5 - $hi) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase - 5 + 0 ));

   
    $p->setcolour("white");
    $p->line( $xbase + $j * $w + $w,      $ysize - ($ybase - 5 - $hi) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase - 5 + 0 ) );
    
    $p->setcolour("red");    
    $p->box({filled => 1}, 
	    $xbase + $j * $w,      $ysize - ($ybase - 5 - $h2) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase - 5 - $h1 ));
    
    my $cc = $c; $cc =~ s/^.+\]//; #print "$c $cc\n";

     $p->setcolour("black");
      $p->setfont("Courrier", 6);
    $p->text( { rotate => 90 }, $xbase + $j * $w + 3*$w/4, $ysize - ($ybase-5 - $hi), $cc);  


    $j ++;
  }

  $p->setlinewidth(0.5);

  $j = 0;
   foreach my $c (@$a_ref_header) {
    
    my $h1 = $hi * ($bins[$j]->[0] - $min_i1 ) / $th ;
    my $h2 = $hi * ($bins[$j]->[1] - $min_i1 ) / $th ;
   
    $p->setcolour("white");
    $p->line( $xbase + $j * $w + $w,      $ysize - ($ybase - 5 - $hi) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase - 5 + 0 ) );
    
    $j ++;
  }

  

  print "$max_i2\t$min_i1\n";

  $p->setfont("Courrier", 8);

  $p->setcolour("black");    
  $p->text( {align => 'center', rotate => 90}, $xbase - 1, $ysize - ($ybase - 5 - $hi + 4), $max_i2);
  $p->text( {align => 'center', rotate => 90}, $xbase - 1, $ysize - ($ybase - 5 - 0   + 1), $min_i1);

}


$p->setcolour("black");
$p->setfont("Courrier", 8);


if ($showstatlabels == 1) {

  my $angle = 60;
  my $off   = 35; my $dec1 = 30; 
  $p->text( {align => 'left', rotate => $angle },  $xbase + $j*$w+35,  $ysize - $ybase, "optimized motif");
  $p->text( {align => 'left', rotate => $angle },  $xbase + $j*$w+85,  $ysize - $ybase, "location");
  $p->text( {align => 'left', rotate => $angle },  $xbase + $j*$w+125, $ysize - $ybase, "MI (bits)");
  $p->text( {align => 'left', rotate => $angle },  $xbase + $j*$w+165, $ysize - $ybase, "z-score");
  $p->text( {align => 'left', rotate => $angle },  $xbase + $j*$w+195, $ysize - $ybase, "robustness");
  #$p->text( {align => 'left', rotate => $angle },  $xbase + $j*$w+225, $ysize - $ybase, "optimization stability");
  $p->text( {align => 'left', rotate => $angle },  $xbase + $j*$w+225, $ysize - $ybase, "position bias");
  $p->text( {align => 'left', rotate => $angle },  $xbase + $j*$w+255, $ysize - $ybase, "orientation bias");
  $p->text( {align => 'left', rotate => $angle },  $xbase + $j*$w+285, $ysize - $ybase, "conservation index"); # +40\
  $p->text( {align => 'left', rotate => $angle },  $xbase + $j*$w+325, $ysize - $ybase, "seed"); # +40
  $p->text( {align => 'left', rotate => $angle },  $xbase + $j*$w+365, $ysize - $ybase, "motif name"); # +40

}

my $min  = undef;
my $max  = undef;
if (defined($densityfile)) {
  $min = 0.0;
  $max = 1.0;
} else {
  $min = -$lp_t_draw;
  $max =  $lp_t_draw;
}

my $cntm = 0;
my $NN   = @$a_ref_header;
my $i    = 0;


foreach my $re (@MOTIFS) {
  
  my $myre = $re;
  if ($STAT{$re}->{RNA} == 1) {
    $myre =~ s/T/U/g;
  }	

  my $r = undef;
  if (!defined($densityfile)) {
    $r = $MATRIX{ $re };  # get the data (row)
  } else {
    $r = $DENSITIES{ $re };
  }

  print "Processing $re ... ";

  #
  #  MOTIF NAME
  #
  my $myii = $MOTIF_NUMBER{ $re };

  print "Outputing motif $myii.eps ... ";

  if ($redoweblogo == 1) {

    my $mo = Sets::myre2wm($myre);
    open OUT, ">$d/$myii.txt" or die "cannot open $d/$myii.txt\n";
    print OUT $mo;
    close OUT;
    
    system("$scriptdir/weblogo/seqlogo -f $d/$myii.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $d/$myii.eps");
    
    if ($dopng == 1) {
      system("$scriptdir/weblogo/seqlogo -f $d/$myii.txt -F PNG  -a -c -M -n -Y -w 5 -h 3 > $d/$myii.png");
    }

  }
  
  my $j = 0;
  foreach my $c (@$r) {
    
    my @col = ();
    if (!defined($colmap)) {
      @col = Sets::interp_general( $c, [0, 0, 204], [255, 255, 0], $min, $max);
    } else {
      @col = Sets::interp_from_matlab_colormap( $c, $A_REF_COLMAP, $min, $max);
    }
    
    $p->setcolour(@col);
    
    $p->box({filled => 1}, 
	    $xbase + $j * $w,      $ysize - ($ybase + $i*$h) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase + ($i*$h+$h)));
    
    $j ++;
  }
  

  #
  # START: integrate Motif LOGO
  #
  my $e  = new PostScript::Simple::EPS(file => "$d/$myii.eps");

  # get height
  my $eh = $e->height;
  my $ew = $e->width;
  # height must be $h, so scale down to $h = k * LO
  $e->scale($h / $eh);
  my $ew_new = int(0.5 + $ew * $h / $eh);

  $p->_add_eps($e, $xbase + $j * $w + 35 - $ew_new / 2,  $ysize - ($ybase + ($i*$h+$h))); 
  
  #
  # END: integrate LOGO
  #
  

  #
  # START : additional info
  #
  $p->setcolour("black");
  $p->setfont("Courrier", 8);


  my $myseed = $STAT{$re}->{SEED};
  if ($STAT{$re}->{RNA} == 1) {
    $myseed =~ s/T/U/g;
  }

  $p->text( {align => 'center' },  $xbase + $j*$w+85,  $ysize - ($ybase + ($i*$h+$h/2)), ($STAT{$re}->{RNA}==1?"3'UTR":"5'"));
  $p->text( {align => 'center' },  $xbase + $j*$w+125, $ysize - ($ybase + ($i*$h+$h/2)), sprintf("%4.3f", $STAT{$re}->{"MI"}));

  $p->text( {align => 'center' },  $xbase + $j*$w+165, $ysize - ($ybase + ($i*$h+$h/2)), sprintf("%3.1f", $STAT{$re}->{"Z"}));
  $p->text( {align => 'center' },  $xbase + $j*$w+195, $ysize - ($ybase + ($i*$h+$h/2)), $STAT{$re}->{R});
#  $p->text( {align => 'center' },  $xbase + $j*$w+225, $ysize - ($ybase + ($i*$h+$h/2)), $STAT{$re}->{S});
  $p->text( {align => 'center' },  $xbase + $j*$w+225, $ysize - ($ybase + ($i*$h+$h/2)), ($STAT{$re}->{DIST}==1?"Y":"-"));
  
  if ($STAT{$re}->{ORIE} == 0) {
    $p->text( {align => 'center' },  $xbase + $j*$w+255, $ysize - ($ybase + ($i*$h+$h/2)), ($STAT{$re}->{ORIE}>0?"Y":"-"));
  } else {
    &draw_arrow($p, $xbase + $j*$w+255, $ysize - ($ybase + ($i*$h+$h/2) - 3), $STAT{$re}->{ORIE});
  }

  $p->text( {align => 'center' },  $xbase + $j*$w+285, $ysize - ($ybase + ($i*$h+$h/2)), (defined($STAT{$re}->{CONS})?$STAT{$re}->{CONS}:'-'));
  my $gaps ;
  my $format=undef;
  if ($myseed =~ /\./)
  {
      $myseed =~ /(\.+)/;
      $gaps=$1;
      $gaps=length($gaps);
      $format = ".($gaps)";
  }
  $myseed =~ s/\.+/$format/ if (defined $format and $gaps>2);
  
  $p->text( {align => 'center' },  $xbase + $j*$w+325, $ysize - ($ybase + ($i*$h+$h/2)), $myseed);

  my $na = $STAT{$re}->{NAME};  # dirty fix
  $na =~ s/^M\d{5}\_//;
  $na =~ s/^J\_M.\d+\_//;
  $na =~ s/\.txt$//;
  $p->text( {align => 'left' },  $xbase + $j*$w+365, $ysize - ($ybase + ($i*$h+$h/2)), $na);

  
  
  #
  # END: additional info
  #

  

  $i ++;
  print "Done.\n";
}


print "Plotting significance boxes.\n";

my $i          = 0;
my $cnt_signif = 0;
my $cntm       = 0;

$p->setlinewidth(0.5);
$p->line( $xbase, $ysize - ($ybase), $xbase + $NN * $w, $ysize - ($ybase));

foreach my $re (@MOTIFS) {
  
  my $r = $MATRIX{$re};
  
  my $j = 0;
  foreach my $c (@$r) {

    my $lp = - Sets::log10( 0.05 / $NN );

    if (abs($c) > $lp) {
      
      my $col = undef;
      if ($c > 0) {
	$col = "red";
      } else {
	$col = "blue";
      }	

      my $dec_top_y = 0;
      my $dec_bot_y = 0;
      
      if (($i != 0) && ($a_ref_clust->[$i-1] ne $a_ref_clust->[$i])) {
	$dec_top_y = 1.5;
      }

      if ($a_ref_clust->[$i] ne $a_ref_clust->[$i+1]) {
	$dec_bot_y = 1.5;
      }


      #print " Motif $re, column $j: logp $c > $lp.\n";
      $p->setcolour( $col );
      $p->box({filled => 0}, 
	      $xbase + $j * $w,      $ysize - ($ybase + $i*$h + $dec_top_y) , 
	      $xbase + $j * $w + $w, $ysize - ($ybase + ($i*$h+$h) - $dec_bot_y));
    }
    $j ++;
  }

  
  $i ++;

}

my $i          = 0;

foreach my $re (@MOTIFS) {
  
  my $r = $MATRIX{$re};

  if ($a_ref_clust->[$i] ne $a_ref_clust->[$i+1]) {
    $p->setcolour( "white" );
    $p->setlinewidth(1.5);    
    $p->line( $xbase, $ysize - ($ybase + $i*$h + $h), $xbase + $NN * $w, $ysize - ($ybase + $i*$h + $h));
    $p->setlinewidth(0.5);
    
  }


  $i ++;
}


$p->setlinewidth(1);



if ($showscale == 1) {
  drawScale(5, $ybase + ($quantized==0?20:0), $min, $max, 50, $p, $xsize, $ysize);
}


if (defined($outeps)) {
  print "Creating $outeps ...";
  $p->output("$outeps");
} else {
  print "Creating $motiffile.eps ...";
  $p->output("$motiffile.eps");
}

print "Done.\n";



if ($ps2pdf == 1) {
  if (defined($outeps)) {
    my $outpdf = $outeps;
    $outpdf =~ s/eps$/pdf/;
    print "Creating PDF $outpdf ... ";
    system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf");
  } else {
    print "Creating PDF $motiffile.pdf ... ";
    system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $motiffile.eps $motiffile.pdf");
  }
  print "Done.\n";
}






sub drawScale {
  my ($x, $y, $min, $max, $res, $p, $xsize, $ysize) = @_;

  my $h = 2; my $w = 20;
  
  my $sep = 0;

  $p->setcolour("black");
  $p->setfont("Courrier", 8);
  $p->text({align => "center"}, $x+$w/2+0, $ysize - ($y - 3), $max);

  $p->setfont("Times", 14);

  if (!defined($densityfile)) {
    $p->text({align => "left", rotate => 90}, $x+$w/2+5, $ysize - ($y - 16), "over-representation");
  }
  
  #print "$min\t$max\n";

  my $t = $max;
  

  for (my $i=0; $i<=$res; $i++) {
    #my @col = interp_general( $t, [255, 255, 0], [204, 0, 0], $min, $max);
    #my @col = Sets::interp_general( $t, [0, 0, 204], [255, 255, 0], $min, $max);
    
    my @col = ();
    if (!defined($colmap)) {
      @col = Sets::interp_general( $t, [0, 0, 204], [255, 255, 0], $min, $max);
    } else {
      @col = Sets::interp_from_matlab_colormap( $t, $A_REF_COLMAP, $min, $max);
    }
    

    #print "$t\t" . join("\t", @col) . "\n";
    
    if ($i >= $res/2) {
      $sep = 10;
    }

    if (defined($densityfile)) {
      $sep = 0;
    }

    $p->setcolour( @col );
    $p->box({filled => 1}, $x, $ysize - ($y + $sep + $i*$h) , $x+$w, $ysize - ($y + $sep + $i*$h + $h));
   
    $t -= ($max - $min) / $res;
    
  }

  $p->setcolour( 'black' );
  $p->setfont("Courrier", 8);

  $p->text({align => "center"}, $x+$w/2+1, $ysize - ($y + $sep + $res*$h + 11), $min);

  

  $p->setfont("Times", 14);
  if (!defined($densityfile)) {
    $p->text({align => "right", rotate => 90}, $x+$w/2+5, $ysize - ($y + $sep + $res*$h + 11 + 10), "under-representation");
  }
  
}





sub interp {
  my ($r, $min, $max) = @_;

  if ($r < $min) {
    $r = $min;
  }

  if ($r > $max) {
    $r = $max;
  }

  #  1 --> red           => "0.8  0    0",
  #  0 --> blue          => "0    0    0.8",
  
  #  0.8 = $max . a + b
  #  0   = $min . a + b

  my $a1 = 0.8 / ($max - $min);
  my $b1 = - $a1 * $min;

  #  0   = $max . a + b
  #  0.8 = $min . a + b

  my $a3 = -0.8 / ($max - $min);
  my $b3 = -$a3 * $max;

  my $c1 = int( 0.5 + 256 * ($a1 * $r + $b1) );
  my $c2 = 0;
  my $c3 = int( 0.5 + 256 * ($a3 * $r + $b3) );
  
  
  return ($c1, $c2, $c3);
  
}


sub sort2DArrayByMax {
  my ($a_ref) = @_;

  my @a = sort sortMax @$a_ref;

  return \@a;
}

sub sortMax {
  #my ($a, $b) = @_;
  
  my $n1 = scalar(@$a);
  my $m1 = -10000;
  my $i1 = undef;
  for (my $i=1; $i<$n1; $i++) {
    if ($a->[$i] > $m1) {
      $m1 = $a->[$i];
      $i1 = $i;
    }    
  }

  my $n2 = scalar(@$b);
  my $m2 = -10000;
  my $i2 = undef;
  for (my $i=1; $i<$n2; $i++) {
    if ($b->[$i] > $m2) {
      $m2 = $b->[$i];
      $i2 = $i;
    }    
  }

  return $i1 <=> $i2;
}




sub binom_test_greater {
  
  my ($n, $k, $p, $p1) = @_;

  my $todo = "$progdir/binom_test_greater $n $k $p 1";
  my $out = `$todo`;  #print "$todo\n";
  $out =~ s/[\n\r]//g;

  my @a = split /\t/, $out, -1;
  
  $out = $a[0];

  $$p1 = $a[1];

  #print "$out\n";
  
  return $out;

}

sub draw_arrow {
  my ($p,$x,$y,$d) = @_;

  
  $p->setcolour( 'black' );

  if ($d == 1) {
    $p->line( $x-5, $y, $x+5, $y);
    $p->line( $x+2, $y+3, $x+5, $y);
    $p->line( $x+2, $y-3, $x+5, $y);
  } else {
    $p->line( $x-5, $y, $x+5, $y);
    $p->line( $x-2, $y+3, $x-5, $y);
    $p->line( $x-2, $y-3, $x-5, $y);
  }

  
}
