#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
my $scriptdir = "$home/PERL_MODULES/SCRIPTS";

use Sets;
use Getopt::Long;
use strict;




my $matrixfile    = undef;
my $logtransform  = 1;
my $qnorm         = 1;
my $varnorm       = 1;
my $cluster       = 1;
my $k             = undef;

if (@ARGV == 0) {
  die "Args: --matrix=FILE --k=INT\n";
}

GetOptions("matrixfile=s"  => \$matrixfile,
	   "k=s"           => \$k);

my $outfile = $matrixfile;
my $infile  = undef;
my $todo    = undef;

if ($logtransform == 1) {

  $infile   = $outfile;
  $outfile .= ".log";

  print "\nLog-transforming the expression values ... ";
  $todo = "perl $scriptdir/expression_log_transform_matrix.pl $infile > $outfile ";
  system($todo);
  print "Done.\n";
}

if ($qnorm == 1) {

  $infile   = $outfile;
  $outfile .= ".qnorm";

  print "\nRunning a quantile normalization using R ... ";

  my $txt = "f <- \"$infile\"
m <- read.csv(f, sep=\"\\t\", row.names=1, header=T)
library(limma)
mn <- normalizeQuantiles(m)
write.table(format(mn, digits=3, trim=T), file=\"$outfile\", sep=\"\\t\", quote=F, col.names=NA, row.names=T)
";
  
  my $tmpfile = Sets::getTempFile("/tmp/Rscript");
  open OUT, ">$tmpfile" or die "Cannot open $tmpfile for writing.\n";
  print OUT $txt;
  close OUT;
  system("R CMD BATCH $tmpfile");

  print "Done.\n";
  
}

if ($varnorm == 1) {

  $infile   = $outfile;
  $outfile .= ".varnorm";

  print "\nVariance normalize matrix ... ";
  $todo = "perl $scriptdir/expression_normalize_rows.pl < $infile > $outfile";
  system($todo);
  print "Done.\n"; 
}


if ($cluster == 1) {

  $infile   = $outfile;
  
  my $nbc = undef;
  if (!defined($k)) {
    $todo = "getnbclusters $infile";
    $nbc = `$todo`; chomp $nbc;
  } else {
    $nbc = $k;
  }

  print "\nCluster into $nbc clusters.\n";
  
  
  $todo = "/usr/local/bin/cluster -g 2 -k $nbc -r 10 -f $infile";
  system($todo) == 0 or die "Could not cluster ... \n";

}
