use lib "$ENV{FIREDIR}/SCRIPTS";
use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

use PostScript::Simple;

# use Data::Dumper;

use Table;
use Sets;

use strict;


my $scriptdir = "$ENV{FIREDIR}/SCRIPTS";
my $progdir   = "$ENV{FIREDIR}/PROGRAMS";

my $colmap    = "$scriptdir/HEATMAPS/cmap_hot.txt";
my $showscale = 1;
my $ps2pdf    = 1;
my $rootdir   = ".";

if (Sets::exist_parameter(\@ARGV, "-colmap") == 1) {
  $colmap   = Sets::get_parameter(\@ARGV, "-colmap");
}
     
my $summaryfile   = Sets::get_parameter(\@ARGV, "-summaryfile");
Sets::defined_exists_or_die($summaryfile, "Please define -summaryfile.\n");

my $matrixfile   = Sets::get_parameter(\@ARGV, "-matrixfile");
Sets::defined_exists_or_die($matrixfile, "Please define -matrixfile.\n");

my $resmatrixfile   = Sets::get_parameter(\@ARGV, "-resmatrixfile");
Sets::defined_exists_or_die($resmatrixfile, "Please define -resmatrixfile.\n");

my $clustfile = undef;
if (Sets::exist_parameter(\@ARGV, "-clustfile") == 1) {
  $clustfile   = Sets::get_parameter(\@ARGV, "-clustfile");
}

my $orderfile = undef;
if (Sets::exist_parameter(\@ARGV, "-orderfile") == 1) {
  $orderfile   = Sets::get_parameter(\@ARGV, "-orderfile");
}


my $motifdir    = Sets::get_parameter(\@ARGV, "-motifdir");
Sets::defined_exists_or_die($motifdir, "Please define -motifdir.\n");

my @MOTIFS = ();

my $ta = Table->new;

#
#  read in the summary file
#
print "Read in summary file ... ";
$ta->loadFile($summaryfile);
my $a_ref_mo = $ta->getArray();
my %STAT      = ();
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
		"CONS"   => $r->[11]);

  $STAT{ $r->[0] } = \%a_tmp;   
  $MOTIF_NUMBER{ $r->[0] } = $cnt ++;
  push @MOTIFS, $r->[0];
}
print "Done.\n";


#
# read in the clusters
#
my $a_ref_clust = undef;
my $h_ref_clust = undef;
if (defined($clustfile)) {
  $ta->loadFile($clustfile);
  $a_ref_clust = $ta->getColumn(1);
  $h_ref_clust = $ta->getIndexKV(0,1);
  my $col         = $ta->getColumn(0);
  @MOTIFS      = @$col;  # over-write
}


if (defined($orderfile)) {
  print "Using orderfile\n";
  $ta->loadFile($orderfile);
  my $col         = $ta->getColumn(0);
  shift @$col if ($col->[0] eq "");
  @MOTIFS      = @$col;  # over-write
}



#
#  read in resmatrix file
#
$ta->loadFile($resmatrixfile);
my $a_ref_resmat = $ta->getArray();

my %RESMATRIX = ();
foreach my $r (@$a_ref_resmat) {
  my $a = shift @$r;
  my $b = shift @$r;
  #print "$a -- $b\n";
  $RESMATRIX{ $a }{ $b } = $r;
  $RESMATRIX{ $b }{ $a } = $r;
  
}


#
#  read in matrix file
#
$ta->loadFile($matrixfile);
my $a_ref_mat = $ta->getArray();

my %MATRIX = ();
foreach my $r (@$a_ref_mat) {
  my $a = shift @$r;
  my $b = shift @$r;
  #print "$a -- $b\n";
  $MATRIX{ $a }{ $b } = $r;
  $MATRIX{ $b }{ $a } = $r;
  
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


print "Now doing the graphical display.\n";

#
#  START DRAWING
#
#

my $xbase = 35;
my $ybase = 115;
my $w     = 30;
my $h     = 30;


my $xsize = $xbase + $w * scalar(@MOTIFS) + 135;
my $ysize = $ybase + $h * scalar(@MOTIFS) + 135; 


print "xsize = $xsize, ysize = $ysize, xbase = $xbase, ybase = $ybase\n";

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



my $min  = -0.01;
my $max  =  0.01;
my $cntm = 0;
my $n    = scalar(@MOTIFS);



#
# add a row of logos
#

my $i = 0;
foreach my $re1 (@MOTIFS) {

  
  #
  # START: integrate Motif LOGO
  #
  
  my $myii = $MOTIF_NUMBER{ $re1 };
  print "$re1\t$myii\n";


  my $d  = $motifdir;  
  my $e  = new PostScript::Simple::EPS(file => "$d/$myii.eps");

  # get height
  my $eh = $e->height;
  my $ew = $e->width;
 
  # height must be $h, so scale down to $h = k * LO
  $e->scale($h / $eh);
  $e->rotate(90);
  my $ew_new = int(0.5 + $ew * $h / $eh);

  # finally add to picture
  $p->_add_eps($e, $xbase + $h+$i*$h + 5,  $ysize - $ybase ); 
  
  #
  # END: integrate LOGO
  #
   
  $p->setcolour("black");
  $p->setfont("Courier-Bold", 10); 
  my $ty = ($STAT{$re1}->{RNA}==0?"5'":"3'UTR");
  $p->text({ rotate => 90, align => 'center'}, $xbase + $h/2+$i*$h+4,  $ysize - ($ybase - $ew_new - 20) , $ty);
  
  $i++;
}
    
$i    = 0;
foreach my $re1 (@MOTIFS) {
  
  my $d       = $motifdir;
  my $j       = 0;
  foreach my $re2 (@MOTIFS) {
    
    my $c = undef;
    
    if ($re1 ne $re2) {
      my $r = $MATRIX{ $re1 }{ $re2 }; 
      if (abs($r->[3]) < 1e-10) {
	$c = 0.0;
      } else {
	$c = $r->[0] * $r->[3] / abs($r->[3]);
      }
    } else {
      $c = 1.0;
    }
      
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
  
  my $myii = $MOTIF_NUMBER{ $re1 };

  
  my $e  = new PostScript::Simple::EPS(file => "$d/$myii.eps");

  # get height
  my $eh = $e->height;
  my $ew = $e->width;
 
  # height must be $h, so scale down to $h = k * LO
  $e->scale($h / $eh);
  my $ew_new = int(0.5 + $ew * $h / $eh);

  # finally add to picture
  $p->_add_eps($e, $xbase + $j * $w + 35 - $ew_new / 2,  $ysize - ($ybase + ($i*$h+$h))); 
  
  #
  # END: integrate LOGO
  #

  $p->setcolour("black");
  #$p->setfont("Courrier", 7); 
  $p->setfont("Courier-Bold", 10); 
  my $ty = ($STAT{$re1}->{RNA}==0?"5'":"3'UTR");
  $p->text({ align => 'center'}, $xbase + $j * $w + $ew_new + 30, $ysize - ($ybase + ($i*$h+$h/2)) , $ty);
  
  
  
  $i ++;

}


if (defined($clustfile)) {

  my $j       = 0;
  my $old_re2 = undef;
  
  $p->setlinewidth(1);
  
  $p->setcolour("black");
  
  $p->line($xbase + $w * 0 ,      $ysize - ($ybase + $j*$h) ,
	   $xbase + $w * $n,      $ysize - ($ybase + $j*$h) );
  
  $p->line($xbase + $w * $j,      $ysize - ($ybase + 0  * $h) ,
	   $xbase + $w * $j,      $ysize - ($ybase + $n * $h) );
  
  
  
  foreach my $re2 (@MOTIFS) {
    
    if ($h_ref_clust->{$re2} != $h_ref_clust->{$old_re2}) {
      
      # fix y, horiz line
      $p->line($xbase + $w * 0 ,      $ysize - ($ybase + $j*$h) ,
	       $xbase + $w * $n,      $ysize - ($ybase + $j*$h) );
      
      # fix x, vertical line
      $p->line($xbase + $w * $j,      $ysize - ($ybase + 0  * $h) ,
	       $xbase + $w * $j,      $ysize - ($ybase + $n * $h) );
      
      
      
    }
    
    $old_re2 = $re2;
    
    $j ++;
  }

  
  
  $p->line($xbase + $w * 0 ,      $ysize - ($ybase + $j*$h) ,
	   $xbase + $w * $n,      $ysize - ($ybase + $j*$h) );
  
  $p->line($xbase + $w * $j,      $ysize - ($ybase + 0  * $h) ,
	   $xbase + $w * $j,      $ysize - ($ybase + $n * $h) );
}


$p->setlinewidth(1);

$i    = 0;
foreach my $re1 (@MOTIFS) {
  my $j     = 0;
  my $rna1  = $STAT{$re1}->{RNA};

  foreach my $re2 (@MOTIFS) {
    my $rna2  = $STAT{$re2}->{RNA};

    #print "$re1 -- $re2\n";

    my $r   = $RESMATRIX{ $re1 }{ $re2 }; 
    my $t   = undef;    

    if (defined($r)) {
      $t = $r->[1];
    } else {
      $t = 10000;
    }

    
    
    
    #print "t = $t\n";
    
    if ($t == 0) {

      if ($rna1 != $rna2) {
	$p->setcolour("green");
      } else {
	if ($rna1 == 0) {
	  $p->setcolour("blue");
	} else {
	  $p->setcolour((255, 102, 255));	  
	}
      }

     
      $p->box({filled => 0}, 
	      $xbase + $j * $w ,      $ysize - ($ybase + $i*$h ) , 
	      $xbase + $j * $w + $w , $ysize - ($ybase + ($i*$h+$h) ));


      # BEWARE: used the shift function twice above !

      my $c1  = 0;
      if (($r->[4] ne "") && ($r->[4] ne "nan") && ($r->[4] <= 100)) { 
	$c1 = 1;
      }
      
      my $c2 = 0;
      if (($r->[6] ne "") && ($r->[6] ne "nan") && ($r->[6] <= 100)) { 
	$c2 = 1;
      }
  
      my $c3 = 0;
      if (($r->[8] ne "") && ($r->[8] ne "nan") && ($r->[8] <= 100)) { 
	$c3 = 1;
      }
      
      my $c4 = 0;
      if (($r->[10] ne "") && ($r->[10] ne "nan") && ($r->[10] <= 100)) { 
	$c4 = 1;
      }

      
      my $c5 = 0;
      if (($r->[12] ne "") && ($r->[12] ne "nan") && ($r->[12] <= 100)) { 
	$c5 = 1;
      }
      
      my $c6 = 0;
      if (($r->[14] ne "") && ($r->[14] ne "nan") && ($r->[14] <= 100)) { 
	$c6 = 1;
      }
      
      $p->setcolour("black");
      $p->setfont("Courrier", 7);

      # noam
  #    if ($c1+$c2+$c3+$c4+$c5+$c6) {      

#	  my $htick = $h/10;
	
#	  if (($c1==1) && ($c2==0)) {					 	      
#	      $p->text({align => "center"}, $xbase + $j * $w + $w/2, $ysize - ($ybase +  $i*$h + 3*$htick    ), "d:i");					 					 
#	  }
#	  if (($c1==0) && ($c2==1)) {					 	      
#	      $p->text({align => "center"}, $xbase + $j * $w + $w/2, $ysize - ($ybase +  $i*$h + 3*$htick    ), "d:u");					 					 
#	  }
#	  if (($c1==1) && ($c2==1)) {					 	      
#	      $p->text({align => "center"}, $xbase + $j * $w + $w/2, $ysize - ($ybase +  $i*$h + 3*$htick    ), "d:i,u");					 					 
#	  }


#	  if (($c3==1) && ($c4==0)) {					 	      
#	      $p->text({align => "center"}, $xbase + $j * $w + $w/2, $ysize - ($ybase +  $i*$h + 6*$htick   ), "o:i");					 					 
#	  }
#	  if (($c3==0) && ($c4==1)) {					 	      
#	      $p->text({align => "center"}, $xbase + $j * $w + $w/2, $ysize - ($ybase +  $i*$h + 6*$htick    ), "o:u");					 					 
#	  }
#	  if (($c3==1) && ($c4==1)) {					 	      
#	      $p->text({align => "center"}, $xbase + $j * $w + $w/2, $ysize - ($ybase +  $i*$h + 6*$htick    ), "o:i,u");					 					 
#	  }


#	  if (($c5==1) && ($c6==0)) {					 	      
#	      $p->text({align => "center"}, $xbase + $j * $w + $w/2, $ysize - ($ybase +  $i*$h + 9*$htick    ), "r:i");					 					 
#	  }
#	  if (($c5==0) && ($c6==1)) {					 	      
#	      $p->text({align => "center"}, $xbase + $j * $w + $w/2, $ysize - ($ybase +  $i*$h + 9*$htick    ), "r:u");					 					 
#	  }
#	  if (($c5==1) && ($c6==1)) {					 	      
#	      $p->text({align => "center"}, $xbase + $j * $w + $w/2, $ysize - ($ybase +  $i*$h + 9*$htick    ), "r:i,u");					 					 
#	  }

#      }

      # OLD
      # $p->text({align => "center"}, $xbase + $j * $w + $w/2, $ysize - ($ybase +  $i*$h + $h/2    ), "$c1$c2$c3$c4$c5$c6");
 
      
      if ($c1 == 1) { # || ($c2 == 1)) {
	
	$p->setcolour("black");
	$p->line(  $xbase + $j * $w + $w/4,
		   $ysize - ($ybase +  $i*$h + $h/2    ) ,
		   $xbase + $j * $w + 3*$w/4,
		   $ysize - ($ybase +  $i*$h + $h/2    ));
	
	$p->line(  $xbase + $j * $w + $w/2,
		   $ysize - ($ybase +  $i*$h + $h/4    ) ,
		   $xbase + $j * $w + $w/2,
		   $ysize - ($ybase +  $i*$h + 3*$h/4    ));
      }
      
      
    }
    $j ++;
    
  }
  $i++;
}


if ($showscale == 1) {
  drawScale(5, $ybase, $min, $max, 50, $p, $xsize, $ysize);
}

print "Creating $matrixfile.eps ...";
$p->output("$matrixfile.eps");
print "Done.\n";



if ($ps2pdf == 1) {
  print "Creating PDF $matrixfile.pdf ... ";
  system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $matrixfile.eps $matrixfile.pdf");
  print "Done.\n";
}






sub drawScale {
  my ($x, $y, $min, $max, $res, $p, $xsize, $ysize) = @_;

  my $h = 2; my $w = 20;
  
  my $sep = 0;

  $p->setcolour("black");
  $p->setfont("Courrier", 8);
  $p->text({align => "center"}, $x+$w/2+0, $ysize - ($y - 3),           $max);
  $p->text({align => "center"}, $x+$w/2+0, $ysize - ($y - 13),           "(Pos)");

  
  #print "$min\t$max\n";

  my $t = $max;
  
  for (my $i=0; $i<$res/2; $i++) {
    my @col = ();
    if (!defined($colmap)) {
      @col = Sets::interp_general( $t, [0, 0, 204], [255, 255, 0], $min, $max);
    } else {
      @col = Sets::interp_from_matlab_colormap( $t, $A_REF_COLMAP, $min, $max);
    }
    $p->setcolour( @col );
    $p->box({filled => 1}, $x, $ysize - ($y + $sep + $i*$h) , $x+$w, $ysize - ($y + $sep + $i*$h + $h));
    $t -= ($max - $min) / $res;
  }

  $p->setcolour("black");
  $p->setfont("Courrier", 8);
  $p->text({align => "center"}, $x+$w/2+0, $ysize - ($y + ($res/2)*$h + 10 ), "0.0");

  $sep = 50;
  
  
  for (my $i=$res-1; $i>=$res/2; $i--) {
    my @col = ();
    if (!defined($colmap)) {
      @col = Sets::interp_general( $t, [0, 0, 204], [255, 255, 0], $min, $max);
    } else {
      @col = Sets::interp_from_matlab_colormap( $t, $A_REF_COLMAP, $min, $max);
    }
    $p->setcolour( @col );
    $p->box({filled => 1}, $x, $ysize - ($y + $sep + $i*$h) , $x+$w, $ysize - ($y + $sep + $i*$h + $h));
    $t -= ($max - $min) / $res;
  }

  $p->setcolour( 'black' );


  
  $p->text({align => "center"}, $x+$w/2+1, $ysize - ($y + $sep + ($res/2-2)*$h-10), "(Neg)");

  $p->text({align => "center"}, $x+$w/2+1, $ysize - ($y + $sep + ($res/2-2)*$h), -$min);
 
  $p->text({align => "center"}, $x+$w/2+0, $ysize - ($y + $sep+$res*$h + 10 ), "0.0");

 
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





