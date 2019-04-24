#!/usr/bin/perl
use lib "$ENV{CHIPSEEQERDIR}";

use Sets;
use strict;

use Getopt::Long;
my $annotation = undef;
my $expfile    = undef;
my $verbose    = 0;
my $add2model  = "";
my $model      = undef;
my $reg        = "log";

if (@ARGV == 0) {
  die "Args: --annotation=FILE --expfile=FILE
where:
--annotation=FILE(S)     annotation files from ChIPseeqerAnnotate
--expfile=FILE           expression file in binary formatx
\n";
}
GetOptions("annotation=s" => \$annotation,
	   "verbose=s"    => \$verbose,
	   "add2model=s"  => \$add2model,
	   "reg=s"        => \$reg,
	   "model=s"      => \$model,
           "expfile=s"    => \$expfile);

open IN, $expfile;
my $l = <IN>; chomp $l;
my @a = split /\t/, $l, -1;
my $explabel = $a[1];
close IN;

my @annofiles  = split /\,/, $annotation;
my @datalabels = ();
my @annofilenames = ();

foreach my $af (@annofiles) {
  open IN, $af or die "Cannot open $af\n";
  $l = <IN>; chomp $l;
  @a = split /\t/, $l, -1;
  shift @a;
  my @labels = @a;
  map s/\-/\./g, @labels;
  push @datalabels, @labels;
  close IN;
  push @annofilenames, Sets::filename($af);
}



# synchronize
my $f      = Sets::filename($expfile);
my $outfile = join(".", @annofilenames) . ".$f";
my $todo = "perl $ENV{CHIPSEEQERDIR}/SCRIPTS/expression_concatenate_matrices.pl " . join(" ", @annofiles) . " $expfile > $outfile";
system($todo) == 0 or die "Cannot exec $todo\n";

if (!defined($model)) {
  $model  =  join(" + ", @datalabels) . " $add2model";
}

my $rscript = "
m <- read.csv(\"$outfile\", header=T, row.names=1, sep=\"\\t\", check.names=T)
attach(data.frame(m))\n";

if ($reg eq "log") {
  $rscript .= "fit <- glm($explabel ~ $model,  family=binomial(link=\"logit\"))\n";
} elsif ($reg eq "linear") {
  $rscript .= "fit <- lm($explabel ~ $model)\n";
}

$rscript .= "print(summary(fit))\n";

if ($verbose == 1) {
  print "$rscript\n";
}
print "$rscript\n";
open OUT, "| R --slave";
print OUT $rscript;
close OUT;





