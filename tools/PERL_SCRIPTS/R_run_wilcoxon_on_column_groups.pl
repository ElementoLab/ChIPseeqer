#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;
use strict;

if (@ARGV == 0) {
  die "args: file numcol1 numcol2\n";
}

my $fi = $ARGV[0];

my $x1 = 1;
my $x2 = $ARGV[1];
my $y1 = $x2+1;
my $y2 = $ARGV[1] + $ARGV[2];

my $tt = "two.sided";
if ($ARGV[3]) {
  $tt = $ARGV[3];
}

my $fo = "$fi.wilcoxpv.$ARGV[1]$tt$ARGV[2]";

my $txt = "
m <- read.csv(\"$fi\", sep=\"\\t\", row.names=1, header=T)
pv <- as.matrix(apply(m, 1, function(x) { wilcox.test(x[$x1:$x2], x[$y1:$y2], alternative=\"$tt\" )\$p.value } ))
write.table(pv, file=\"$fo\", sep=\"\\t\", row.names=T, col.names=NA, quote=F)
";

my $ft = Sets::getTempFile("/tmp/Rscript");
Sets::writeText($txt, $ft);

system("R CMD BATCH $ft") == 0 or die "Cannot exec R script ?\n";

system("cat $fo");

unlink $ft;


