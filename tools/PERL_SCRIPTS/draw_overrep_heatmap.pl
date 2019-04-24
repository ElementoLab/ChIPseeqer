#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use strict;
use Matrix2EPS;
use AggloClust;

use Getopt::Long;
use Table;

my $matrix          = undef;
my $clustcols       = 0;
my $clustrows       = 0;
my $h               = 10;
my $w               = 10;
my $coldesc         = undef;
my $rowdesc         = undef;
my $reordercols     = undef;
my $suffix          = undef;
my $minmax          = 3;
my $here            = 0;
my $saverownames    = undef;
my $savecolnames    = undef;
my $draw            = undef;
my $outfile         = undef;
my $xbase           = 80;
my $ybase           = 150;
my $sortcolsbygene  = undef;
my $cmap            = "/Users/olivier/PROGRAMS/MYTREEVIEW/HEATMAPS/cmap2_tab.txt";
my $algoclust       = undef;
my $xright          = 250;
my $headerangle     = undef;
my $drawhorizscale  = undef;
my $headerfontsize  = undef;
my $rownamefontsize = undef;
my $savematrix      = undef;
my $euclidean       = 1;
my $drawsignifboxes = 0;
my $drawonlysignif  = 0;
my $tsignif         = 0.01;
my $rowdrawmotifs   = undef;
my $coldrawmotifs   = undef;
my $distance        = undef;
my $drawheader      = undef;
my $maxtext         = "Enrichment";
my $mintext         = "Depletion";
my $scaleres        = 120;
my $xscale          = undef;
my $sortrowsbymax   = undef;
my $sortrowsbymin   = undef;

if (@ARGV == 0) {
  die "args: --matrix --clustcols --clustrows

Other options:
--drawonlysignif=INT
--headerfontsize=INT
--savecolnames=FILE
\n";
}

GetOptions("matrix=s"          => \$matrix,
	   "clustcols=s"       => \$clustcols,
	   "clustrows=s"       => \$clustrows,
	   "h=s"               => \$h,
	   "w=s"               => \$w,
	   "draw=s"            => \$draw,
	   "drawheader=s"      => \$drawheader,
	   "reordercols=s"     => \$reordercols,
	   "cmap=s"            => \$cmap,
	   "sortcolsbygene=s"  => \$sortcolsbygene,
	   "rowdesc=s"         => \$rowdesc,
	   "suffix=s"          => \$suffix,
	   "minmax=s"          => \$minmax,
	   "here=s"            => \$here,
	   "xbase=s"           => \$xbase,
	   "ybase=s"           => \$ybase,
	   "xright=s"          => \$xright,
	   "outfile=s"         => \$outfile,
	   "algoclust=s"       => \$algoclust,
	   "saverownames=s"    => \$saverownames,
	   "savecolnames=s"    => \$savecolnames,
	   "savematrix=s"      => \$savematrix,
	   "coldesc=s"         => \$coldesc,
	   "headerangle=s"     => \$headerangle,
	   "drawhorizscale=s"  => \$drawhorizscale,
	   "headerfontsize=s"  => \$headerfontsize,
	   "rownamefontsize=s" => \$rownamefontsize,
	   "euclidean=s"       => \$euclidean,
	   "drawsignifboxes=s" => \$drawsignifboxes,
	   "drawonlysignif=s"  => \$drawonlysignif,
	   "tsignif=s"         => \$tsignif,
	   "sortrowsbymax=s"   => \$sortrowsbymax,
	   "sortrowsbymin=s"   => \$sortrowsbymin,	   
	   "rowdrawmotifs=s"   => \$rowdrawmotifs,
           "coldrawmotifs=s"   => \$coldrawmotifs,
	   "distance=s"        => \$distance,
	   "maxtext=s"         => \$maxtext,
	   "mintext=s"         => \$mintext,
	   "scaleres=s"        => \$scaleres,
	   "xscale=s"          => \$xscale);


my $ma = Matrix2EPS->new;

if (!defined($outfile)) {
  $outfile = "$matrix";
  if (defined($suffix)) {
    $outfile .= ".$suffix";
  }
  $outfile .= ".eps";
}


if ($here == 1) {
  $outfile = Sets::filename($outfile);
}


$ma->setOutputFileName($outfile);
my $mymin = -$minmax;
my $mymax = $minmax;
if ($minmax =~ /\,/) {
  my @a = split /\,/, $minmax;
  $mymin = $a[0];
  $mymax = $a[1];
}


$ma->colnameAngle(90);

# EMP

#$ma->drawHeaderArrows("High in NB, low in GC", "Relative expression", "Low in GC, high in GC"); 
#$ma->drawHeader(0);
$ma->setMinAbsLP(3);
$ma->drawOnlyPos(1); # only show enrichment
$ma->setMinText($mintext);
$ma->setScaleRes($scaleres);

$ma->drawBoxUp(-Sets::log10($tsignif));
$ma->drawSignifBoxes($drawsignifboxes);
$ma->drawOnlySignif($drawonlysignif);
$ma->setFont("Arial");
$ma->setMin($mymin);
$ma->setMax($mymax);
$ma->setMaxText($maxtext);
$ma->setVerbose(1);
$ma->setH($h);
$ma->setW($w);
$ma->setAlgoClust($algoclust) if (defined($algoclust));
#$ma->addColumnDesc($ARGV[1]);
$ma->setColMap($cmap);

if (defined($headerangle)) {
  $ma->colnameAngle($headerangle);
}

$ma->drawRowNames("right");

#$ma->setNeg(1);
$ma->loadMatrix($matrix, 1, 1);

# $ma->setMatrix(\@M, \@motifs, \@regions);

if (defined($coldesc)) {
  print "Adding column info from $coldesc\n";
  $ma->addColumnDesc($coldesc);
}

if (defined($rowdesc)) {
  print "Adding row info from $rowdesc\n";
  $ma->addRowDesc($rowdesc);
}

if ($clustrows == 1) {
  $ma->clusterRows($distance);
} elsif (defined($sortrowsbymax)) {
  $ma->sortRowsByMax($sortrowsbymax);
} elsif (defined($sortrowsbymin)) {
  $ma->sortRowsByMin($sortrowsbymin);
}

if ($clustcols == 1) {
  $ma->clusterColumns($distance);
  #$ma->clusterColumns();
}

if ($reordercols == 1) {
  print "Reordering cols .. \n";
  $ma->reorderColumnsUsingDesc();
}

if (defined($sortcolsbygene)) {
  $ma->sortColumnsByGene($sortcolsbygene);
}

$ma->setXbase($xbase);
$ma->setYbase($ybase);
$ma->setXright($xright);

#$ma->loadIclustPartition($ARGV[1]);
#$ma->reorderUsingPartition(1,1);

$ma->regexpColNames('s/\.txt$//');
$ma->regexpColNames('s/.calls.ORF//');
$ma->regexpColNames('s/motif_//');
$ma->regexpColNames('s/.txt.ORF.noamb//');
$ma->regexpColNames('s/.txt.+$//');
$ma->regexpColNames('s/genelist_//');

$ma->regexpRowNames('s/.txt.+$//');
$ma->regexpRowNames('s/.txt$//');
$ma->regexpRowNames('s/genelist_//');
$ma->regexpRowNames('s/^M\d+_//');
$ma->regexpRowNames('s/^J_MA\d+_//');


if (defined($rowdrawmotifs) && ($rowdrawmotifs == 1)) {
  $ma->rowDrawMotifs(1);
}



if (defined($coldrawmotifs) && ($coldrawmotifs == 1)) {
  $ma->colDrawMotifs(1);
}

$ma->setRownameFontSize((defined($rownamefontsize)?$rownamefontsize:$h+1));
$ma->setHeaderFontSize((defined($headerfontsize) ?$headerfontsize :$w+1));  # was $h+7 for Ari

#$ma->loadColumnsPartition($ARGV[1]);


$ma->draw();
#$ma->addColumnClusters();

if ($drawhorizscale == 0) {
  $ma->drawScale(20,250);
} else {
  # $xbase + @{$a_ref_M->[0]} * $w / 2 - 100 - 100, $ybase + @$a_ref_M*$h+ $h*8/4
  my $dim = $ma->getMatrixDimensions();

  print "$dim->[0]\t$dim->[1]\n";

  $ma->setScaleRes($scaleres); # tmp
  #$ma->drawHorizontalScale($xbase + $dim->[1] * $w / 2 - 100 - 100,
  #			   $ybase + $dim->[0] * $h + $h*2);
  
#  $ma->drawHorizontalScale($xbase + $w/2 ,
#  			   $ybase + $dim->[0] * $h + 20 + 20);
  
  my $x = $xbase + $dim->[1] * $w / 2 - 2 * $scaleres / 2 + 50;
  if (defined($xscale)) {
    $x = $xscale;
  }
  $ma->drawHorizontalScale($x,
  			   $ybase + $dim->[0] * $h + 20 + 20);

  
}

$ma->output();
$ma->pdfify();



if (defined($saverownames)) {
  $ma->saveOrderedRowNames($saverownames);
}

if (defined($savecolnames)) {
  $ma->saveOrderedColNames($savecolnames);
}

if (defined($savematrix)) {
  $ma->saveMatrix($savematrix);
}


$outfile =~ s/\.eps/\.pdf/;

print "Created $outfile\n";

if (defined($draw)) {

  my $todo = "$draw $outfile";
  system($todo);
}
