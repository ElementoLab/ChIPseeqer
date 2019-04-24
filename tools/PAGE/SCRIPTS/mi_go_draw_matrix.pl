my $pagedir ;
BEGIN{
    if ((!$ENV{PAGEDIR}) || ($ENV{PAGEDIR} eq '')) {
	$pagedir="./" ;
	print "The PAGEDIR environment variable is not set. It is set to default.\n";
    }
    else{
	$pagedir = $ENV{PAGEDIR};
    }
}

use lib "$pagedir/SCRIPTS";
use lib "$pagedir/SCRIPTS/PostScript-Simple-0.07/lib";

my $programdir = $pagedir."/PROGRAMS" ;
my $scriptdir  = $pagedir."/SCRIPTS" ;

use Table;
use Sets;
use Getopt::Long;
use PostScript::Simple;
use AggloClust;
use strict;
use Data::Dumper ;

my $expdata             = undef;
my $pvaluematrixfile    = undef;
my $expfile             = undef;
my $datafile            = undef;
my $colmap              = "$scriptdir/HEATMAPS/cmap_dens.txt";
my $cluster             = undef;
my $sortrowsbyphase     = 1;
my $max_p               = undef;
my $minmax_lp           = 3;
my $quantized           = 1;
my $min                 = undef;
my $max                 = undef;
my $xsize               = undef;
my $xscale              = 10 ;
my $yscale              = 50 ;
my $scalefont           = undef ;
my $h                   = undef ;
my $w                   = undef ;
my $draw_sample_heatmap = "false" ;
my $order               = 0 ;
my $fontsize            = 20;
my $arrowcentertext     = undef;
my $arrowlefttext       = undef;
my $arrowrighttext      = undef;
my $onecat              = undef;
my $font                = "Arial";
my $drawframes          = 1;
my $catstoignore        = undef;
my $discretize          = 1;
my $collabels           = undef;
my $ybase               = 100;
my $onlysignif          = 0;
my $contheader          = 1;
my $horizscale          = 0;
my $mirror              = 0;
my $pmax                = 0.01;
my $ignore_non_signif   = 0;
my $rotateheader        = 90;

if (@ARGV == 0) {
 die "Usage: perl mi_go_draw_matrix.pl  --pvaluematrixfile=FILE --expfile=FILE --max_p=P\n";
}

GetOptions ('pvaluematrixfile=s'    => \$pvaluematrixfile,
            'cluster=s'             => \$cluster,
	    'expfile=s'             => \$expfile,
	    'datafile=s'            => \$datafile,
	    'quantized=s'           => \$quantized,
	    "expdata=s"             => \$expdata,
	    'onecat=s'              => \$onecat,
	    'colmap=s'              => \$colmap,
	    "rotateheader=s"        => \$rotateheader,
	    'minmax_lp=s'           => \$minmax_lp,
	    'min=s'                 => \$min,
	    'max=s'                 => \$max,
	    'font=s'                => \$font,
	    'xsize=s'               => \$xsize,
	    'xscale=s'              => \$xscale,
	    'yscale=s'              => \$yscale,
	    'ybase=s'               => \$ybase,
	    'onlysignif=s'          => \$onlysignif,
	    'scalefont=s'           => \$scalefont,
	    'collabels=s'           => \$collabels,
	    'arrowcentertext=s'     => \$arrowcentertext,
	    'arrowlefttext=s'       => \$arrowlefttext,
	    'arrowrighttext=s'      => \$arrowrighttext,
	    'drawframes=s'          => \$drawframes,
	    'h=s'                   => \$h,
	    'w=s'                   => \$w,
	    'max_p=s'               => \$max_p,
	    'fontsize=s'            => \$fontsize,
	    'catstoignore=s'        => \$catstoignore,
	    'draw_sample_heatmap=s' => \$draw_sample_heatmap,
	    'order=s'               => \$order,
	    'contheader=s'          => \$contheader,
	    'horizscale=s'          => \$horizscale,
	    'mirror=s'              => \$mirror,
	    'ignore_non_signif=s'   => \$ignore_non_signif,
	    "pmax=s"                => \$pmax);

#
# creating the summary file
#
my $file = substr($expfile, rindex($expfile, '/')) ;
my $dir = $expfile."_PAGE/" ;
my $ta = Table->new;

print "Reading matrix ... ";

#
#  read in the matrix file
#
$ta->loadFile($pvaluematrixfile);

# get an 2D array
my $a_ref_M      = $ta->getArray();

if ($mirror == 1) {

  my @NEWM = ();
  
  foreach my $r (@$a_ref_M) {
    
    my @row = ();
    $row[0] = $r->[0];
    for (my $i=@$r-1, my $j=1; $i>=1; $i--, $j++) {
      $row[$j] = $r->[$i];
    }
    push @NEWM, \@row;

  }

  @$a_ref_M = @NEWM;
  
}

# header
my $a_ref_H      = shift @$a_ref_M; shift @$a_ref_H;
if (!defined($max_p)) {
  $max_p = 0.05 / @$a_ref_H;
}
print "Done.\n";

if (defined($onecat)) {
  shift @$a_ref_M;
}

if (defined($catstoignore)) {
  
  print "Ignoring cats ... \n";
  
  my $h_ref_i = Sets::getIndex($catstoignore);

  my @NEWM = ();
  foreach my $r (@$a_ref_M) {
    my @a = split /\ /, $r->[0];
    print "$a[0]\n";
    if (defined($h_ref_i->{$a[0]})) {
      print "Ignoring $a[0]\n";
    } else {
      push @NEWM, $r;
    }
  }
  
  $a_ref_M = \@NEWM;
}

if (defined($ignore_non_signif) && ($ignore_non_signif == 1)) {

  #print "I am here\n"; <STDIN>;
  my @NEWM = ();
  foreach my $r (@$a_ref_M) {  
    # go through entries in current row

    my $hassignif = 0;
    for (my $j=1; $j<@$r; $j++) {

      # get log pvalues for over and under rep
      my ($lpo,$lpu) = $r->[$j] =~ /^(.+?)\/(.+)$/;
      # get a single value, pos if over-rep, neg if under-rep
      my $lpt = Sets::log10($pmax); #/@$r);
      #print "$lp > $lpt?\n";
      if (abs($lpo) > abs($lpt)) {
	$hassignif = 1;
      }
      
    } # loop over cols    
    if ($hassignif == 1) {
      print join(" ", @$r) . "\n"; 
      push @NEWM, $r;
    }
  } # loop over rows
  
  $a_ref_M = \@NEWM;
  

}

if (defined($cluster) && (@$a_ref_M > 2)) {

  print "Cluster rows .. ";
  $cluster = Sets::min(scalar(@$a_ref_M), $cluster);

  my $ac = AggloClust->new;
  
  my @dist = ();
  my $n    = @$a_ref_M;
  for (my $i=0; $i<$n-1; $i++) {
    $dist[$i][$i] = 0;

    for (my $j=$i+1; $j<$n; $j++) {
      
      my @a1 = @{ $a_ref_M->[$i] }; shift @a1;
      
      for (my $k=0; $k<@a1; $k++) {
	
	my ($lpo,$lpu) = $a1[$k] =~ /^(.+?)\/(.+)$/;	
	if (abs($lpo) > abs($lpu)) {
	  $a1[$k] = -$lpo;
	} else {
	  $a1[$k] = $lpu;
	}
	
	if ($discretize == 1) {
	  if ($a1[$k] > -Sets::log10(0.01)) {
	    $a1[$k] = 1;
	  } elsif ($a1[$k] < Sets::log10(0.01)) {
	    $a1[$k] = -1;
	  } else {
	    $a1[$k] = 0;
	  }
	} 
       
      }
    
      my @a2 = @{ $a_ref_M->[$j] }; shift @a2;
      
      #print join("\t", @a2) . "\n";
      
      for (my $k=0; $k<@a2; $k++) {
	my ($lpo,$lpu) = $a2[$k] =~ /^(.+?)\/(.+)$/;	
	if (abs($lpo) > abs($lpu)) {
	  $a2[$k] = -$lpo;
	} else {
	  $a2[$k] = $lpu;
	}
	if ($discretize == 1) {
	  if ($a2[$k] > -Sets::log10(0.01)) {
	    $a2[$k] = 1;
	  } elsif ($a2[$k] < Sets::log10(0.01)) {
	    $a2[$k] = -1;
	  } else {
	    $a2[$k] = 0;
	  }
	}
      }      
      #print join("\t", @a2) . "\n";
    
      
      
      $dist[$i][$j] = 1 - Sets::pearson(\@a1, \@a2);
      $dist[$j][$i] = $dist[$i][$j];
    }
  }
  
  
  $ac->setDistanceMatrix(\@dist);
 $ac->setUseCorr(0);

 #$ac->setMaxNbClusters($cluster);
 #my $a_ref_c = $ac->agglomerate_using_avg_linkage();
 
 $ac->agglomerate_using_max_linkage();
 my $a_ref_o = $ac->getDFSOrder();

 my @NEWMAT = ();
 foreach my $c (@$a_ref_o) {
   #print join(" ", @$c); print "\n";
   #foreach my $i (@$c) {
   push @NEWMAT, $a_ref_M->[$c];
   #}
 }
 $a_ref_M = \@NEWMAT;

 print "Done.";
}

#
# load color map
#
my $A_REF_COLMAP = undef;
if (defined($colmap)) {

 $ta->setDelim("[\t\ ]");
 $ta->loadFile($colmap);
 $A_REF_COLMAP = $ta->getArray();
 $ta->setDelim('\t');
}

# font size
if (!defined($scalefont)) {
  $scalefont = $fontsize;
}



#
#  START DRAWING
#
#
my %CENTROIDS ;
my @samples ;
my @clusters ;
my $nbsamples ;
my $minimum ;
my $maximum ;
my %ARRAYS ;

#
# left and top margins
#
my $xbase        = 60;

$ybase           = 100 if ($quantized==0) ;

# height of each entry in the matrix
$h            = 30 if ! (defined $h) ;

# size of the EPS image to be generated
my $ysize        = $ybase + $h * scalar(@$a_ref_M) + 100 + 200;  #Sets::max(600,$ybase + $h * scalar(@$a_ref_M) + 100);
my $h_h = 10 ;
if ($quantized == 1 and defined $datafile and $draw_sample_heatmap eq "true")
{
    my $todo = "perl $scriptdir/draw_clusters_and_gocats.pl --expfile=$expfile --data=$datafile --pvaluematrixfile=$pvaluematrixfile" ;
    system("$todo") == 0 or die "system failed: $?";
    
    open CL, "< $dir/sample_centroids.txt" or die;
    print "$dir/sample_centroids.txt" ;
    
    chomp ;
    
    $minimum = 1000;
    $maximum = 0;
    my $l = <CL> ;
    chomp ($l) ;
    @samples = split(/\t/, $l) ;
    shift (@samples) ;
    while(<CL>)
    {
	chomp ;
	my ($cluster, @a) = split(/\t/, $_) ;
	$nbsamples = $#a+1 ;
	for (my $i=0; $i<$nbsamples; $i++) 
	{
	    $CENTROIDS{$cluster}{$samples[$i]} = $a[$i] ;
	    push(@{$ARRAYS{$samples[$i]}}, $a[$i]) ;
	}
	my $m1 = Sets::minInArray(\@a) ;
	my $m2 = Sets::maxInArray(\@a) ;
	$minimum = $m1 if ($m1<$minimum) ;
	$maximum = $m2 if ($m2>$maximum) ;
	push (@clusters, $cluster) ;
    }
    $ysize += $h_h * $nbsamples ;

    for (my $j=0; $j<$nbsamples; $j++){
	my $average = Sets::average($ARRAYS{$samples[$j]}) ;
	my $std = Sets::stddev($ARRAYS{$samples[$j]}) ;
	for (my $i=0 ; $i<=$#clusters ; $i++){
	    $CENTROIDS{$clusters[$i]}{$samples[$j]} = ($CENTROIDS{$clusters[$i]}{$samples[$j]}-$average)/$std;
	}
    }

    my @CEN ;
    my @NEWMAT = ();
    my @HEADER = () ;
    my @enr ;
    for (my $i=0 ; $i<@$a_ref_H ; $i++){
	my $c = $a_ref_H->[$i] ;
	my $delta = $CENTROIDS{$c}{$samples[1]}-$CENTROIDS{$c}{$samples[0]} ;
	$enr[$i]->{delta} = $delta ;
	$enr[$i]->{index} = $i+1 ;
    }
    @enr = sort {$b->{delta} <=> $a->{delta}} (@enr) ;

    for (my $i=0 ; $i<@$a_ref_M ; $i++){
	$NEWMAT[$i][0] = $a_ref_M->[$i][0] ;
    }
    my $cnt = 1 ;
    foreach my $e (@enr){
        my $j = $e->{index} ;
	print $j, "\t" ;
	for (my $i=0 ; $i<@$a_ref_M ; $i++){
	    $NEWMAT[$i][$cnt] = $a_ref_M->[$i][$j] ;
	}
	$HEADER[$cnt-1] = $a_ref_H->[$j-1] ;
	$cnt++ ;
    }
    $a_ref_M = \@NEWMAT ;
    $a_ref_H = \@HEADER ;
}

if ($order ==1){
    my @NEWMAT = ();
    my @enr ;
    for (my $i=0 ; $i<@$a_ref_M ; $i++){
	my $r = $a_ref_M->[$i] ;
	my $min = 0 ;
	my $pos = 1 ;
	for (my $j=1 ; $j<@$r ; $j++){
	    my ($lpo,$lpu) = $r->[$j] =~ /^(.+?)\/(.+)$/ ;
	    if ($lpo<$min){
		$min = $lpo ;
		$pos = $j ;
	    }
	}
	$enr[$i]->{ind} = $i ;
	$enr[$i]->{pos} = $pos ;
    }
    @enr = sort {$b->{pos} <=> $a->{pos}} (@enr) ;
    foreach my $e (@enr){
	my $i = $e->{ind} ;
	push @NEWMAT, $a_ref_M->[$i];
    }
    $a_ref_M = \@NEWMAT ;
}


$xsize           = 1400 if ! (defined $xsize);

# width of each entry in the matrix
if (! defined($w)) {
  my $canvasize = Sets::min(500,$xsize/2.0);
  $w = int( 0.5 + ((1 * $canvasize) + 5)  / @$a_ref_H )-1;
}
print "Start drawing\n";

my $p = new PostScript::Simple(xsize     => $xsize,
			       ysize     => $ysize,
                              colour    => 1,
                              eps       => 1,
                              units     => "pt");

$p->setlinewidth(0.5);

my $f = Sets::filename($expfile) ;
#
# show header (cluster indices)
#
if ($quantized == 1) {

    if (defined $datafile and $draw_sample_heatmap eq "true") {
	my @color = ();
	#my $cmap = "$scriptdir/HEATMAPS/c.txt" ;
	my $cmap = undef ;
	my $REF_COLMAP = undef;
	for (my $j=0; $j<$nbsamples; $j++) 
	{
	    for (my $i=0 ; $i<@$a_ref_H ; $i++)
	    {
		if (defined($cmap)) 
		{
		    $ta->setDelim(" ");
		    $ta->loadFile($cmap);
		    $REF_COLMAP = $ta->getArray();
		    $ta->setDelim('\t');
		}
		my @color = () ;
		if ($CENTROIDS{$a_ref_H->[$i]}{$samples[$j]}>0){
		    @color = Sets::interp_general( $CENTROIDS{$a_ref_H->[$i]}{$samples[$j]},[255, 0, 0], [255, 255, 0], 0, 2);
		}
		else{
		    @color = Sets::interp_general( $CENTROIDS{$a_ref_H->[$i]}{$samples[$j]},[0, 0, 0], [255, 0, 0], -2, 0)
		}
		#my @color = Sets::interp_general( $CENTROIDS{$a_ref_H->[$i]}{$samples[$j]},[0, 255, 0], [255, 0, 0], -2, 2);
		#my @color = Sets::interp_from_matlab_colormap( $CENTROIDS{$a_ref_H->[$i]}{$samples[$j]},$REF_COLMAP, -2, 2);
		
		$p->setcolour(@color);
		$p->box({filled => 1},
			$xbase + $i * $w,      $ysize - ($ybase + $j*$h_h) ,
			$xbase + $i * $w + $w, $ysize - ($ybase + ($j*$h_h+$h_h)));
		
	    }
	    $p->setfont("Times", 8);
	    $p->setcolour("black");
	    $p->text({align => "left", rotate => 0}, $xbase + scalar(@$a_ref_H) * $w + 10, $ysize - ($ybase + $j*$h_h+$h_h/2+4), $samples[$j]);
	    print $samples[$j], "\n" ;
	}

	drawScale($xscale+20, $ybase , -2, 2, 50, $p, $xsize, $ysize, $cmap, $REF_COLMAP, "", "", 8, 0.3, 20);
	$ybase = $ybase - ($nbsamples*$h_h) + $h_h*6 + 10;
      } # if datafile provided


    #
    # added by OE to draw expression heatmap / looks similar to what Hani coded above 
    #
    if (defined($expdata)) {  
      
      #
      # load heatmap
      #
      %CENTROIDS = ();
      open IN, $expdata or die "Cannot open $expdata\n";
      my $l = <IN>; chomp $l;
      my @clids = split /\t/, $l;
      shift @clids;
      my @samples = ();
      while (my $l = <IN>) {
	chomp $l;	
	my @a = split /\t/, $l, -1;
	my $n = shift @a;
	push @samples, $n;
	for (my $i=0; $i<@a; $i++) {
	  $CENTROIDS{ $clids[$i] }{ $n } = $a[$i]
	}
      }
      close IN;
      $nbsamples = @samples;
      
      # need expression cmap
      $ta->setDelim("[\t\ ]");
      my $cmap = "$ENV{HEATMAPDIR}/colormaps/rbg_cmap.csv";
      if (! -e $cmap) {
	die "Cannot find $cmap\n";
      }
      $ta->loadFile("$cmap");
      my $REF_COLMAP = $ta->getArray();
      $ta->setDelim('\t');
      
      my $expdataminmax = 1;

      my @color = ();
      $h_h = $h;
      for (my $j=0; $j<$nbsamples; $j++) {
	for (my $i=0 ; $i<@$a_ref_H ; $i++) {
	  
	  my $c = $CENTROIDS{$a_ref_H->[$i]}{$samples[$j]};
	  @color = Sets::interp_from_matlab_colormap( $c, $REF_COLMAP, -$expdataminmax, $expdataminmax);
	  
	  $p->setcolour(@color);
	  $p->box({filled => 1},
		  $xbase + $i * $w,      $ysize - ($ybase + $j*$h_h) ,
		  $xbase + $i * $w + $w, $ysize - ($ybase + ($j*$h_h+$h_h)));
	  
	}
	$p->setfont("Times", $h_h);
	$p->setcolour("black");
	$p->text({align => "left", rotate => 0}, $xbase + scalar(@$a_ref_H) * $w + 10, $ysize - ($ybase + $j*$h_h+$h_h/2+10), $samples[$j]);
	print $samples[$j], "\n" ;
      }
      
      drawScale($xscale+20, $ybase , -$expdataminmax, $expdataminmax, 50, $p, $xsize, $ysize, 1, $REF_COLMAP, "", "", 8, 0.3, 20);
      $ybase = $ybase + ($nbsamples*$h_h) + $h_h; #+ $h_h*6 + 10;
      
    } # oe addition





    my $h_ref_lab = undef;
    if (defined($collabels)) {
      
      my $ta1 = Table->new;
      $ta1->loadFile($collabels);
      $h_ref_lab = $ta1->getIndexKV(0,1);
      
    }

    

    # drawing header

    $p->setcolour("black");
    $p->setfont($font, $fontsize);
    
    for (my $j=0; $j<@$a_ref_H; $j++) {

      my $cc = $a_ref_H->[$j]; 
      
      if (defined($h_ref_lab) && defined($h_ref_lab->{$cc})) {
	$cc = $h_ref_lab->{$cc};
      }
      
      $p->text( { rotate => $rotateheader }, $xbase + $j * $w + 3 * $w/4, $ysize - ($ybase-3), $cc);
    }
  
# end (if quantized == 1)  
} else {
    


  if (defined($arrowcentertext)) {
    
    $p->setcolour("black");
    
    my $nc = @$a_ref_H;
    
    my $x1 = $xbase;
    my $x2 = $x1 + $nc * $w;

    my $y1 = $ysize - ($ybase - $h * 1.5 - 10);
    my $y2 = $y1;

    $p->line($x1, $y1, $x2, $y2);

    # arrow left
    $p->line($x1, $y1, $x1+5, $y1-5);
    $p->line($x1, $y1, $x1+5, $y1+5);

    # arrow right
    $p->line($x2, $y2, $x2-5, $y2-5);
    $p->line($x2, $y2, $x2-5, $y2+5);


    $p->setfont($font, 10);
    $p->text({ align => "left"}, $x1+5, $y1+2, $arrowlefttext);
    $p->text({ align => "right" }, $x2-5, $y2+2, $arrowrighttext);

    $p->text({ align => "centre" }, ($x2+$x1)/2, $y2+20, $arrowcentertext);

  }


  if ($contheader == 1) {
    
    # graphical header for continuous data
    
    my $min_i1 =  1000000;
    my $max_i2 = -1000000;
    my @bins   = ();
    foreach my $c (@$a_ref_H) {
      my ($i1, $i2) = $c =~ /\[(.+?)\ (.+?)\]/;    
      $min_i1 = $i1 if ($i1 < $min_i1);
      $max_i2 = $i2 if ($i2 > $max_i2);
      my @a_tmp = ($i1, $i2); push @bins, \@a_tmp;
    }
    
    my $th = $max_i2 - $min_i1;
    my $hi = $h * 1.5;
    
    my $j = 0;
    foreach my $c (@$a_ref_H) {
      
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
      
      
      $j ++;
    }
    
    $p->setlinewidth(0.5);
    
    $j = 0;
    foreach my $c (@$a_ref_H) {
      
      my $h1 = $hi * ($bins[$j]->[0] - $min_i1 ) / $th ;
      my $h2 = $hi * ($bins[$j]->[1] - $min_i1 ) / $th ;
      
      $p->setcolour("white");
      $p->line( $xbase + $j * $w + $w,      $ysize - ($ybase - 5 - $hi) , 
		$xbase + $j * $w + $w, $ysize - ($ybase - 5 + 0 ) );
      
      $j ++;
    }
    
    
    
    #print "$max_i2\t$min_i1\n";
    
    $p->setfont("Courrier", 8);
    
    $p->setcolour("black");    
    $p->text( {align => 'left'}, $xbase + 1 + scalar(@$a_ref_H) * $w , $ysize - ($ybase - 5 - $hi + 4), $max_i2);
    $p->text( {align => 'right'}, $xbase - 1, $ysize - ($ybase - 5 - 0   + 1), $min_i1);
    
  
  }

}


$p->setcolour("black");
$p->setfont("Courrier", 8);


#
# draw (i,j) p-value matrix itself
#

# 
my $outlabelfile = "$expfile"."\_PAGE/$f".".summary.labels.txt";
open OUTLABELS, ">$outlabelfile";


# set min max p-values
if ($minmax_lp =~ /\,/) {
  my @a = split /\,/, $minmax_lp;
  $min  = $a[0];
  $max  = $a[1];
} else {
  $max =  $minmax_lp if !(defined $max);
  $min = -$minmax_lp if !(defined $min);
}

my @col = ();
my @go ;
for (my $i=0; $i<@$a_ref_M; $i++) {
  # get row
  my $r  = $a_ref_M->[$i];
  
  # get GO description
  my $go = shift @$r;
  print OUTLABELS "$go\n";
  
  # 
  $go =~ s/^(.+?)\ //;
  $go = "$go, $1" if ($go ne $1) ;
  
  my @goc = split //, $go;
  if (($goc[1] ne uc($goc[1])) && ($goc[2] ne uc($goc[2]))) {
   $goc[0] = uc($goc[0]);
 }
  
  $go = join("", @goc);
  
  print "$go\n";
  push(@go, $go) ;
  
  # go through entries in current row
 for (my $j=0; $j<@$r; $j++) {

   # get log pvalues for over and under rep
   my ($lpo,$lpu) = $r->[$j] =~ /^(.+?)\/(.+)$/;

   # get a single value, pos if over-rep, neg if under-rep
   my $lp = undef;
   if (abs($lpo) > abs($lpu)) {
     $lp = -$lpo;
   } else {
     $lp = $lpu;
   }
   my $lpt = Sets::log10($pmax); #/@$r);

   # create appropriate color
   my @col = ();  #$colmap = undef;


   if (($onlysignif == 1) && (abs($lp) < abs($lpt))) {
     @col = "black";

   } elsif (($onlysignif == 2) && ( (abs($lp) < abs($lpt)) || ($lp < 0)  )  ) {
     @col = "black";

   } else {
     
     if (!defined($colmap)) {
       if ($lp<0)
	 {
	   @col = Sets::interp_general( $lp, [0, 255, 0], [0, 0, 0], $min, 0);
	 }
       else
	 {
	   @col = Sets::interp_general( $lp, [0, 0, 0], [255, 0, 0], 0, $max);
	 }
     } else {
       @col = Sets::interp_from_matlab_colormap( $lp, $A_REF_COLMAP, $min, $max);
     }
     
   }
   
   # draw the matrix entry with appropriate color
   $p->setcolour(@col);
   $p->box({filled => 1},
           $xbase + $j * $w,      $ysize - ($ybase + $i*$h) ,
           $xbase + $j * $w + $w, $ysize - ($ybase + ($i*$h+$h)));

 }

 $p->setfont($font, $fontsize);
 $p->setcolour("black");
 $p->text({align => "left", rotate => 0}, $xbase + @$r * $w + 10, $ysize - ($ybase + $i*$h+3*$h/4), $go);


}

close OUTLABELS;



if ($drawframes == 1) {

  for (my $i=0; $i<@$a_ref_M; $i++) {

    

    # get row
    my $r  = $a_ref_M->[$i];

    ## get GO description
    #my $go = shift @$r;

    # go through entries in current row
    for (my $j=0; $j<@$r; $j++) {
      
      # get log pvalues for over and under rep
      my ($lpo,$lpu) = $r->[$j] =~ /^(.+?)\/(.+)$/;
      
      # get a single value, pos if over-rep, neg if under-rep
      my $lp = undef;
      if (abs($lpo) > abs($lpu)) {
	$lp = -$lpo;
      } else {
	$lp = $lpu;
      }

      my $lpt = Sets::log10($pmax); #/@$r);
      
      if (($onlysignif == 2) && ( (abs($lp) < abs($lpt)) || ($lp < 0)  )  ) {
	next;
      }


      if (abs($lp) > abs($lpt)) {
	
	
	# create appropriate color
	my @col = ();  #$colmap = undef;
	if (!defined($colmap)) {

	  if ($lp<0) {
	    @col = Sets::interp_general( $min, [0, 255, 0], [0, 0, 0], $min, 0);
	  } else {
	    @col = Sets::interp_general( $max, [0, 0, 0], [255, 0, 0], 0, $max);
	  }
	
	} else {
	  @col = Sets::interp_from_matlab_colormap( ($lp<0?$min:$max), $A_REF_COLMAP, $min, $max);
	}
	
	# draw the matrix entry with appropriate color
	$p->setcolour(@col);
	$p->box({filled => 0},
		$xbase + $j * $w,      $ysize - ($ybase + $i*$h) ,
		$xbase + $j * $w + $w, $ysize - ($ybase + ($i*$h+$h)));
	
      } # end if significant

    }

  }
  
}





#
# draw scale bar
#
if ($horizscale == 0) {
  drawScale($xscale, $ybase + $yscale , $min, $max, 50, $p, $xsize, $ysize, $colmap, $A_REF_COLMAP, "Enrichment", "Depletion", $scalefont, 2, 20);
} else {
  drawHorizontalScale($xbase + @{$a_ref_M->[0]} * $w / 2 - 100 - 100, $ybase + @$a_ref_M*$h+ $h*8/4 , $min, $max, 100, $p, $xsize, $ysize, $colmap, $A_REF_COLMAP, "Enrichment", "Depletion", $scalefont, 20, 2, 1, undef);
}

my $outeps = "$expfile"."\_PAGE/$f".".summary.eps" ;
my $outpdf = "$expfile"."\_PAGE/$f".".summary.pdf" ;
my $ps2pdf = 1;


# output EPS file
print "Outputing EPS file $outeps\n";
$p->output("$outeps");

# convert to PDF
print "Convert to PDF $outpdf\n";
system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf");

print "Finished.\n";
exit(0);


sub drawScale {
 my ($x, $y, $min, $max, $res, $p, $xsize, $ysize, $colmap, $A_REF_COLMAP, $upper_text, $lower_text, $scalefont, $h, $w) = @_;

 my $sep = 0;

 $p->setcolour("black");
 $p->setfont($font, $scalefont);
 $p->text({align => "center"}, $x+$w/2+0, $ysize - ($y - 3), $max);

 $p->setfont($font, $scalefont);

 $p->text({align => "left", rotate => 90}, $x+$w/2+$scalefont/3, $ysize - ($y - 3 - $scalefont - 3), $upper_text);


 my $t = $max;


 my $fpv = 0;
 for (my $i=0; $i<=$res; $i++) {

   my $dec = ($max - $min) / $res;

   #my $iii = int($t);
   #if (($iii >= 2) && (abs($t-$iii) < $dec)) {
   #  $p->setcolour( 'black' );
   #  my $smallfont = 8;
   #  $p->setfont($font, $smallfont);
   #  $p->text({align => "left"}, $x+$w+2, 
   #	      $ysize - ($y + $sep + $i*$h+$smallfont/2),
   #	      "p<1e-$iii");
   #  $fpv ++;
   #}
   

   my @col = () ;
   if (!defined($colmap)) {
       if ($i>$res/2)
       {
	   @col = Sets::interp_general( $t, [0, 0, 0], [255, 0, 0], $min, 0);
       }
       else
       {
	   @col = Sets::interp_general( $t, [255, 0, 0], [255, 255, 0], 0, $max);
	 }
       
   } else {

     @col = Sets::interp_from_matlab_colormap( $t, $A_REF_COLMAP, $min, $max);
   }

   $p->setcolour( @col );
   $p->box({filled => 1}, $x, $ysize - ($y + $sep + $i*$h) , $x+$w, $ysize - ($y + $sep + $i*$h + $h));
   $t -= ($max - $min) / $res;
 }

 $p->setcolour( 'black' );
 $p->setfont($font, $scalefont);
 $p->text({align => "center"}, $x+$w/2+1, $ysize - ($y + $sep + $res*$h + 1 * $scalefont), $min);
 $p->setfont($font, $scalefont);
 $p->text({align => "right", rotate => 90}, $x+$w/2+$scalefont/3, $ysize - ($y + $sep + $res*$h + 1.25 * $scalefont), $lower_text);


}




sub drawHorizontalScale {
 my ($x, $y, $min, $max, $res, $p, $xsize, $ysize, $colmap, $A_REF_COLMAP, $upper_text, $lower_text, $scalefont, $h, $w, $drawp, $minabslp) = @_;

 my $sep = 0;

 if ($drawp == 0) {
   # print MAX
   $p->setcolour("black");
   $p->setfont($font, $scalefont);
   $p->text({align => "right"}, $x-3, $ysize - ($y+$h-3), $max);
   
   # print TEXT that goes with MAX
   $p->setfont($font, $scalefont);
   $p->text({align => "right"}, $x - 3 - $scalefont, $ysize - ($y+$h-3), $upper_text);
 }

 my $t = $max;

 my $fpv = 0;
 for (my $i=0; $i<=$res; $i++) {

   my $dec = ($max - (defined($minabslp)?$minabslp:$min)) / $res;

   if ($drawp == 1) {
     my $iii = int($t);
     if (($iii >= 0) && (abs($t-$iii) < $dec)) {
       $p->setcolour( 'black' );
       my $smallfont = $scalefont*(3/4);
       $p->setfont($font, $smallfont);
       
       my $txt = undef;
       if ($iii == 0) {
	 $txt = "p=1";	 
       } elsif ($fpv == 0) {
	 $txt = "p<=1e-$iii";
       } else {
	 $txt = "p=1e-$iii";
       }
       

       $p->text({align => "left", rotate => 45}, 
		$x+$i*$w+$smallfont/4,
		$ysize - ($y - 2),
		$txt);
       $fpv ++;
     }
     

   }
   

   my @col = () ;
   if (!defined($colmap)) {
       if ($i>$res/2)
       {
	   @col = Sets::interp_general( $t, [0, 0, 0], [255, 0, 0], $min, 0);
       }
       else
       {
	   @col = Sets::interp_general( $t, [255, 0, 0], [255, 255, 0], 0, $max);
       }
   } else {
     @col = Sets::interp_from_matlab_colormap( $t, $A_REF_COLMAP, $min, $max);
   }

   $p->setcolour( @col );
   $p->box({filled => 1}, $x + $i*$w, $ysize - ($y) , $x + ($i+1)*$w  , $ysize - ($y + $h));
   $t -= $dec;
 }

 #$p->setcolour( 'black' );
 #$p->setfont($font, $scalefont);
 #$p->text({align => "center"}, $x+$w/2+1, $ysize - ($y + $sep + $res*$h + 1 * $scalefont), $min);

 #$p->setfont($font, $scalefont);
 #$p->text({align => "right", rotate => 90}, $x+$w/2+$scalefont/3, $ysize - ($y + $sep + $res*$h + 1.25 * $scalefont), $lower_text);

 if ($drawp == 0) {
   # print MIN
   $p->setcolour("black");
   $p->setfont($font, $scalefont);
   $p->text({align => "left"}, $x+3+$res*$w, $ysize - ($y+$h-3), $min);
   
   # print TEXT that goes with MIN
   $p->setfont($font, $scalefont);
   $p->text({align => "left"}, $x + 3+$res*$w + $scalefont, $ysize - ($y+$h-3), $lower_text);
 } else {

   $p->setcolour("black");
   $p->setfont($font, $scalefont);
   $p->text({align => "center"}, $x+$res*$w/2, $ysize - ($y+$h+$scalefont), "Enrichment significance");

   if (defined($minabslp)) {
     
     $p->setcolour("black");
     my $xx = $x + $res*$w + $h*2;
     $p->box({filled => 1}, $xx, $ysize - ($y) , $xx+$h, $ysize - ($y + $h));

     $p->setcolour("black");
     $p->setfont($font, $scalefont);

     $p->text($xx+$h+3, $ysize - ($y + $h - 3), "Non-significant"); 

   }

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
