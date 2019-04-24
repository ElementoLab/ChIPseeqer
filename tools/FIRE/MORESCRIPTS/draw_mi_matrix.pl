#
#  input summary, mimatrix
#
use lib "$ENV{FIREDIR}/SCRIPTS";
use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

use PostScript::Simple;
use AggloClust;

use Data::Dumper;

use Table;
use Sets;

use strict;

my $colmap    = "$ENV{FIREDIR}/HEATMAPS/cmap_hot.txt";
my $showscale = 1;
my $ps2pdf    = 1;
my $rootdir   = ".";
my $posonly   = 0;


if (Sets::exist_parameter(\@ARGV, "-posonly") == 1) {
  $posonly   = Sets::get_parameter(\@ARGV, "-posonly");
}

if (Sets::exist_parameter(\@ARGV, "-colmap") == 1) {
  $colmap   = Sets::get_parameter(\@ARGV, "-colmap");
}


my $clust = undef;
if (Sets::exist_parameter(\@ARGV, "-clust") == 1) {
  $clust = Sets::get_parameter(\@ARGV, "-clust");
}

my $clustcols = undef;
if (Sets::exist_parameter(\@ARGV, "-clustcols") == 1) {
  $clustcols = Sets::get_parameter(\@ARGV, "-clustcols");
}


my $sortmotifsbyphase = undef;
if (Sets::exist_parameter(\@ARGV, "-sortmotifsbyphase") == 1) {
  $sortmotifsbyphase = Sets::get_parameter(\@ARGV, "-sortmotifsbyphase");
}

my $w     = 30;
if (Sets::exist_parameter(\@ARGV, "-w") == 1) {
  $w   = Sets::get_parameter(\@ARGV, "-w");
}

my $matrixfile   = Sets::get_parameter(\@ARGV, "-matrixfile");
Sets::defined_exists_or_die($matrixfile, "Please define -matrixfile.\n");

my $motiftypefile = undef;
if (Sets::exist_parameter(\@ARGV, "-motiftypefile")) {
  $motiftypefile = Sets::get_parameter(\@ARGV, "-motiftypefile");
}

my $consfile = undef;
if (Sets::exist_parameter(\@ARGV, "-consfile")) {
  $consfile = Sets::get_parameter(\@ARGV, "-consfile");
}

my $restrictfile = undef;
if (Sets::exist_parameter(\@ARGV, "-restrictfile")) {
  $restrictfile = Sets::get_parameter(\@ARGV, "-restrictfile");
}




my $d           = "TMP";
mkdir $d if (! -e $d);

my $ta          = Table->new;


#
# read motif types
#
my %TYPES = ();
if (defined($motiftypefile)) {
  
  $ta->loadFile($motiftypefile);
  my $a_ref_ty = $ta->getArray();
  foreach my $r (@$a_ref_ty) {
    $TYPES{ $r->[0] } = $r->[1];
  }

}


my $h_ref_cons = undef;
if (defined($consfile)) {
  $ta->loadFile($consfile);
  $h_ref_cons = $ta->getIndexKV(0,1);
}

#
#  read in motif restrict file
#
my $h_ref_restrict = undef;
if (defined($restrictfile)) {
  $h_ref_restrict = Sets::getIndex($restrictfile);
}


#
#  read in the matrix file
#
$ta->loadFile($matrixfile);
my $a_ref_mat = $ta->getArray();
if (defined($h_ref_restrict)) {
  my @a = ();
  push @a, shift(@$a_ref_mat);
  push @a, shift(@$a_ref_mat);

  foreach my $r (@$a_ref_mat) {
    push @a, $r if (defined($h_ref_restrict->{$r->[0]}));
    print "retain $r->[0]\n";
  }
  $a_ref_mat = \@a;
  
  print "Got " . scalar(@$a_ref_mat) . " motifs\n";
}

my @MATRIX  = ();
my @MOTIFS  = ();
my @NAMES   = ();
my @ATYPES  = ();

my $a_ref_h = shift @$a_ref_mat;  # header (file names)
my $a_ref_c = shift @$a_ref_mat;  # experiment type (can be file name)
shift @$a_ref_h; # motif regexp 
shift @$a_ref_c; # motif name


#
#  reorganize rows if necessary
#

if (defined($sortmotifsbyphase)) {
  
  my @dist = ();
  my $n    = @$a_ref_mat;
  for (my $i=0; $i<$n; $i++) {
    my @a1 = @{ $a_ref_mat->[$i] }; my $re = shift @a1; shift @a1; 
    my $j = Sets::indexMaxInArray(\@a1);
    push @dist, $j;
    #print "$re\t$j\n";
  }

  my $a_ref_c = Sets::order(\@dist);

  my @NEWMAT = ();
  foreach my $c (@$a_ref_c) {
    push @NEWMAT, $a_ref_mat->[$c];
  }
  $a_ref_mat = \@NEWMAT;

}

#
# cluster the rows
#
if (defined($clust)) {
  
  my $ac = AggloClust->new;

  my @dist = ();
  my $n    = @$a_ref_mat;
  for (my $i=0; $i<$n-1; $i++) {
    $dist[$i][$i] = 0;
    for (my $j=$i+1; $j<$n; $j++) {
      my @a1 = @{ $a_ref_mat->[$i] }; shift @a1; shift @a1; 
      my @a2 = @{ $a_ref_mat->[$j] }; shift @a2; shift @a2; 
      $dist[$i][$j] = 1 - Sets::pearson(\@a1, \@a2);      
      $dist[$j][$i] = $dist[$i][$j]; 
    }
  }

  $ac->setDistanceMatrix(\@dist);
  $ac->setMaxNbClusters($clust);
  my $a_ref_c = $ac->agglomerate_using_avg_linkage();

  my @NEWMAT = ();
  foreach my $c (@$a_ref_c) {
    print join(" ", @$c); print "\n";
    foreach my $i (@$c) {
      push @NEWMAT, $a_ref_mat->[$i];
    }
  }
  $a_ref_mat = \@NEWMAT;
  
}


#
# cluster the cols
#
if (defined($clustcols)) {

  my $ac = AggloClust->new;

  my $a_ref_mat_t = Sets::transpose($a_ref_mat);
  my $l1 = shift @$a_ref_mat_t; # take off regexp
  my $l2 = shift @$a_ref_mat_t; # take off motif name

  my @dist = ();
  my $n    = @$a_ref_mat_t;
  for (my $i=0; $i<$n-1; $i++) { 
    $dist[$i][$i] = 0;
    for (my $j=$i+1; $j<$n; $j++) {
      my @a1 = @{ $a_ref_mat->[$i] }; #shift @a1; #shift @a1; 
      my @a2 = @{ $a_ref_mat->[$j] }; #shift @a2; #shift @a2; 
      $dist[$i][$j] = 1 - Sets::pearson(\@a1, \@a2);      
      $dist[$j][$i] = $dist[$i][$j]; 
    }
  }

  $ac->setDistanceMatrix(\@dist);
  $ac->setMaxNbClusters($clustcols);
  my $a_ref_c = $ac->agglomerate_using_avg_linkage();

  my @NEWMAT = ();
  my @NEWL1  = ();
  my @NEWL2  = ();

  foreach my $c (@$a_ref_c) {
    # print join(" ", @$c); print "\n";
    foreach my $i (@$c) {
      push @NEWMAT, $a_ref_mat_t->[$i];
      push @NEWL1 , $l1->[$i];
      push @NEWL2 , $l2->[$i];
    }
  }

  @NEWMAT = (@NEWL1, @NEWL2, @NEWMAT);

  $a_ref_mat = Sets::transpose(\@NEWMAT);

}


foreach my $r (@$a_ref_mat) {

  my $a = shift @$r;
  my $t = shift @$r;
  push @MOTIFS, $a;
  push @MATRIX, $r;
  push @NAMES,  $t;
  if (defined($motiftypefile)) {
    push @ATYPES, $TYPES{$a};
  } else {
    push @ATYPES, 0;
  }	
}

my $cnt_names = 0;
foreach my $r (@NAMES) {
  $cnt_names ++ if ($r ne '-');
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


my $ybase = 120;
foreach my $s (@$a_ref_h) {
  my $l = length($s) * 6;
  if ($l > $ybase) {
    $ybase = $l;
  }
}

my $h     = 30;
my $xsize = $xbase + $w * scalar(@{$MATRIX[0]}) + 250;
my $ysize = $ybase + $h * scalar(@MOTIFS)       + 135; 


print "xsize = $xsize, ysize = $ysize, xbase = $xbase, ybase = $ybase\n";

my $p = new PostScript::Simple(xsize     => $xsize,
			       ysize     => $ysize,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");

$p->setlinewidth(0.5);



my $min  =  100000;
my $max  =  -100000;
my $cntm = 0;
my $k    = scalar(@MOTIFS);
my $n    = scalar(@{$MATRIX[0]});


for (my $i=0; $i<$k; $i++) {
  for (my $j=0; $j<$n; $j++) {
    $max = $MATRIX[$i][$k] if ($MATRIX[$i][$k] > $max);
    $min = $MATRIX[$i][$k] if ($MATRIX[$i][$k] < $min);
  }
}
    

my $min  =  -10;
my $max  =   10;

#
# add a row of logos
#




my $i       = 0;
my $cnt     = 0;
foreach my $re (@MOTIFS) {
  
  for (my $j=0; $j<$n; $j++) {
    
    my $c = $MATRIX[$i][$j];
    #my $c = rand(0.1);
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
  }
  
  
  #
  # make motif logo
  #
  my $mo = Sets::myre2wm($re);
  if ($ATYPES[$i] == 1) {
    $mo =~ s/T/U/g;
  }
  open OUT, ">$d/$cnt.txt" or die "cannot open $d/$cnt.txt\n";
  print OUT $mo;
  close OUT;
  print "Outputing $re\n";
  system("$ENV{FIREDIR}/SCRIPTS/weblogo/seqlogo -f $d/$cnt.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $d/$cnt.eps");

  #
  # START: integrate Motif LOGO
  #  
  my $e  = new PostScript::Simple::EPS(file => "$d/$cnt.eps");

  # get height
  my $eh = $e->height;
  my $ew = $e->width;
 
  # height must be $h, so scale down to $h = k * LO
  $e->scale($h / $eh);
  #$e->rotate(90);
  my $ew_new = int(0.5 + $ew * $h / $eh);

  # finally add to picture
  $p->_add_eps($e, $xbase + $n*$w + 5,  $ysize - ( $ybase + $i*$h+$h)); 
  
  #
  # END: integrate LOGO
  #
  
  $p->setcolour("black");
  $p->setfont("Courrier", 8);
  
  $p->text({ align => 'center'}, $xbase + $n*$w+$ew_new+25, $ysize - ( $ybase + $i*$h+$h/2), ($ATYPES[$i]==0?"5'":"3'UTR"));  

  my $offset = 0;
  if (defined($consfile)) {
    $p->text({ align => 'left'}, $xbase + $n*$w+$ew_new+50, $ysize - ( $ybase + $i*$h+$h/2), $h_ref_cons->{$re});    
    $offset += 35;
  }

  if ($cnt_names > 0) {
    my $tmp = $NAMES[$i]; $tmp =~ s/\.txt//; $tmp =~ s/^.+\_//g;
    $p->text({ align => 'left'}, $xbase + $n*$w+$ew_new+$offset+50, $ysize - ( $ybase + $i*$h+$h/2), $tmp);
  }

  $i++; $cnt++;
}
    


$p->setcolour("black");
$p->setfont("Courrier", 8);
my $j_last = 0;
push @$a_ref_h, "dummy";

#$p->setcolour((255,255,0));    
$p->line($xbase + 0 * $w, $ysize - $ybase, $xbase + (@$a_ref_h-1) * $w, $ysize - $ybase);


#$p->setcolour((255,255,0));    
$p->line($xbase + 0 * $w, $ysize - ($ybase + $k * $h), $xbase + (@$a_ref_h-1) * $w, $ysize - ($ybase + $k * $h));


#$p->setcolour((255,255,0));    
$p->line($xbase + 0 * $w, $ysize - $ybase, $xbase + 0 * $w, $ysize - ($ybase + $k * $h));

$p->setfont("Courrier", 10);

for (my $j=1; $j<scalar(@$a_ref_h); $j++) {

  if ($a_ref_h->[$j-1] ne $a_ref_h->[$j]) {
    $p->setcolour("black");    

    my $txt = $a_ref_h->[$j-1];
    $txt =~ s/\_/\ /g; $txt =~ s/\.txt//; $txt =~ s/\.sgn//;
    $p->text({ align => 'left', rotate => 90 }, $xbase + 4 + $j_last * $w + (($j - $j_last) * $w) / 2, $ysize - ($ybase  - 8), $txt);
    #$p->setcolour((255,255,0));    
    #$p->line($xbase + $j * $w, $ysize - $ybase, $xbase + $j * $w, $ysize - ($ybase + $k * $h));
    $j_last = $j;
  }
}


 
$p->setcolour("black");
$p->setfont("Courrier", 8);

$p->text({ align => 'left', rotate => 60 }, $xbase + $n*$w+35, $ysize - ($ybase - 8), "Optimized motif");  

$p->text({ align => 'left', rotate => 60 }, $xbase + $n*$w+75, $ysize - ($ybase - 8), "Location");  

my $offset = 0;
if (defined($consfile)) {
  $p->text({ align => 'left', rotate => 60 }, $xbase + $n*$w+110, $ysize - ($ybase - 8), "Conservation index");  
  $offset += 30;
}

if ($cnt_names > 0) {
  $p->text({ align => 'left', rotate => 60 }, $xbase + $n*$w+$offset+110, $ysize - ($ybase - 8), "Motif name");
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

  if ($posonly == 0) {
    $p->text({align => "left", rotate => 90}, $x+$w/2+3, $ysize - ($y - 15), "Positive correlation");
  }
  
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


  #$p->text({align => "center"}, $x+$w/2+0, $ysize - ($y + ($res/2)*$h + 10 ), '0');

  $sep = 0;
  
  if ($posonly == 0) {    
    for (my $i=$res/2; $i<$res; $i++) {
      my @col = ();
      if (!defined($colmap)) {
	@col = Sets::interp_general( $t, [0, 0, 204], [255, 255, 0], $min, $max);
      }  else {
	@col = Sets::interp_from_matlab_colormap( $t, $A_REF_COLMAP, $min, $max);
      } 
      $p->setcolour( @col );
      $p->box({filled => 1}, $x, $ysize - ($y + $sep + $i*$h) , $x+$w, $ysize - ($y + $sep + $i*$h + $h));
      $t -= ($max - $min) / $res;
    }     
    $p->setcolour( 'black' );           
    #$p->text({align => "center"}, $x+$w/2+1, $ysize - ($y + $sep + ($res/2-2)*$h-10), "(Neg)");    
    #$p->text({align => "center"}, $x+$w/2+1, $ysize - ($y + $sep + ($res/2-2)*$h), -);    
    $p->text({align => "center"}, $x+$w/2+0, $ysize - ($y + $sep+$res*$h + 10 ), "-10");
    $p->text({align => "right", rotate => 90}, $x+$w/2+3, $ysize - ($y + $sep + $res * $h + 15), "Negative correlation");    
    
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

  my $todo = "$rootdir/binom_test_greater $n $k $p 1";
  my $out = `$todo`;  #print "$todo\n";
  $out =~ s/[\n\r]//g;

  my @a = split /\t/, $out, -1;
  
  $out = $a[0];

  $$p1 = $a[1];

  #print "$out\n";
  
  return $out;

}

