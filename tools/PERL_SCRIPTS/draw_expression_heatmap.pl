#!/usr/bin/perl

use lib "$ENV{HOME}/PERL_MODULES";

use strict;
use Matrix2EPS;
use AggloClust;

use Getopt::Long;
use Table;

if ($ENV{HEATMAPDIR} eq "") {
  die "You must define \$HEATMAPDIR, using export HEATMAPDIR=/path/to/dir\n";
}

my $matrix           = undef;
my $clustcols        = 0;
my $clustrows        = 0;
my $h                = 10;
my $w                = 10;
my $coldesc          = undef;
my $rowdesc          = undef;
my $reordercols      = undef;
my $suffix           = undef;
my $minmax           = 3;
my $here             = 0;
my $saverownames     = undef;
my $savecolnames = undef;
my $draw             = undef;
my $outfile          = undef;
my $xbase            = 80;
my $ybase            = 100;
my $xright           = undef;
my $sortcolsbygene   = undef;
my $cmap             = "$ENV{HEATMAPDIR}/colormaps/rbg_cmap.csv";
my $sortrowsbycorrel = undef;
my $sortrowsbymax    = undef;
my $sortrowsbymin    = undef;
my $drawrownames     = "right";
my $algoclust        = undef;
my $mintext          = "Down";
my $maxtext          = "Up";
my $distance         = undef;
my $colfile          = undef;
my $yscale           = 150;
my $drawscale        = 1;
my $sortrowsusingfile= undef;
my $righttable       = undef;
my $changenamedcolsign    = undef;
my $neg              = 0;

if (@ARGV == 0) {
  die "args: --matrix=FILE --clustcols=INT --clustrows==INT --h=INT --w=INT
Addditional options:
--saverownames=FILE      to save sorted row names to a file
--algoclust=STR          min, max or avg
--neg=INT                1 to change sign of matrix

\n";
}

GetOptions("matrix=s"            => \$matrix,
	   "clustcols=s"         => \$clustcols,
	   "clustrows=s"         => \$clustrows,
	   "h=s"                 => \$h,
	   "w=s"                 => \$w,
	   "righttable=s"        => \$righttable,
	   "draw=s"              => \$draw,
	   "reordercols=s"       => \$reordercols,
	   "savecolnames=s"      => \$savecolnames,
	   "cmap=s"              => \$cmap,
	   "sortcolsbygene=s"    => \$sortcolsbygene,
	   "sortrowsbycorrel=s"  => \$sortrowsbycorrel,
	   "sortrowsusingfile=s" => \$sortrowsusingfile,
	   "sortrowsbymax=s"     => \$sortrowsbymax,
	   "sortrowsbymin=s"     => \$sortrowsbymin,
	   "rowdesc=s"           => \$rowdesc,
	   "suffix=s"            => \$suffix,
	   "drawrownames=s"      => \$drawrownames,
	   "minmax=s"            => \$minmax,
	   "maxtext=s"           => \$maxtext,
	   "mintext=s"           => \$mintext,
	   "here=s"              => \$here,
	   "xbase=s"             => \$xbase,
	   "xright=s"            => \$xright,
	   "algoclust=s"         => \$algoclust,	
	   "ybase=s"             => \$ybase,
	   "outfile=s"           => \$outfile,
	   "saverownames=s"      => \$saverownames,
	   "coldesc=s"           => \$coldesc,
	   "colfile=s"           => \$colfile,
	   "distance=s"          => \$distance,
	   "neg=s"               => \$neg,
	   "drawscale=s"         => \$drawscale,
	   "changenamedcolsign=s"=> \$changenamedcolsign,
	   "yscale=s"            => \$yscale);


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
$ma->setMin(-$minmax);
$ma->setMax($minmax);
$ma->setMinText($mintext);
$ma->setMaxText($maxtext);
$ma->setVerbose(1);
$ma->setH($h);
$ma->setW($w);
$ma->setAlgoClust($algoclust) if (defined($algoclust));

#$ma->addColumnDesc($ARGV[1]);
$ma->setColMap($cmap);

$ma->drawFrames(0);

$ma->drawRowNames($drawrownames);


if (defined( $changenamedcolsign)) {
  $ma->changeNamedColSign( $changenamedcolsign);
}
if ($neg == 1) {
 $ma->setNeg(1);
}
$ma->loadMatrix($matrix, 1, 1);

if (defined($righttable)) {
  $ma->loadRightTable($righttable);
}


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
} elsif (defined($sortrowsbycorrel)) {
  print "Sorting rows by how they correlate with $sortrowsbycorrel\n";
  $ma->sortRowsUsingCorrelationWithGene($sortrowsbycorrel);
} elsif (defined($sortrowsbymax)) {
  $ma->sortRowsByMax($sortrowsbymax);
} elsif (defined($sortrowsbymin)) {
  $ma->sortRowsByMin($sortrowsbymin);
} elsif (defined($sortrowsusingfile)) {
  $ma->sortRowsUsingFile($sortrowsusingfile);
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

if (defined($colfile)) {
  $ma->reorderColsUsingColList($colfile);
}

$ma->setXbase($xbase);
$ma->setYbase($ybase);
#$ma->setYbase(250);

if (defined($xright)) {
  $ma->setXright($xright);
}

#$ma->loadIclustPartition($ARGV[1]);
#$ma->reorderUsingPartition(1,1);

$ma->regexpRowNames('s/\.txt$//');
$ma->regexpColNames('s/.activity//');
$ma->regexpColNames('s/^.*point\ //');
$ma->regexpColNames('s/Biological Replicate\ //g');

#$ma->colDrawMotifs(1);

$ma->setRownameFontSize($h);
$ma->setHeaderFontSize($w);

#$ma->loadColumnsPartition($ARGV[1]);


$ma->draw();
#$ma->addColumnClusters();
if ($drawscale == 1) {
  $ma->drawScale(20,$yscale);
}

print "Creating heatmap in $outfile (PostScript format) ... ";
$ma->output();
print "Done\n";

$outfile =~ s/\.eps/\.pdf/;

print "Making PDF ";
$ma->pdfify();
if (-e $outfile) {
  print " $outfile ... Done\n";
} else {
  print "Error. You probably do not have ps2pdf installed on your machine.\n";
}


if (defined($saverownames)) {
  $ma->saveOrderedRowNames($saverownames);
}


if (defined($savecolnames)) {
  $ma->saveOrderedColNames($savecolnames);
}

if (defined($draw)) {
  $outfile =~ s/\.eps/\.pdf/;
  my $todo = "$draw $outfile";
  system($todo);
}
