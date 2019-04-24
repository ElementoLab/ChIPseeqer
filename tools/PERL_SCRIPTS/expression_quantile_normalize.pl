#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

my $infile   = $ARGV[0];
my $outfile .= "$infile.qnorm";

print "Running a quantile normalization using R ... ";

my $txt = "f <- \"$infile\"
m <- read.csv(f, sep=\"\\t\", row.names=1, header=T, check.names = FALSE)
";

if ($ARGV[1] ne "") {
  $txt .= "m <- log(m+1)\n";
}

$txt .= "
library(limma)
mn <- normalizeQuantiles(m)
write.table(mn, file=\"$outfile\", sep=\"\\t\", quote=F, col.names=NA, row.names=T)
";

my $tmpfile = Sets::getTempFile("/tmp/Rscript");
open OUT, ">$tmpfile" or die "Cannot open $tmpfile for writing.\n";
print OUT $txt;
close OUT;
system("R CMD BATCH $tmpfile");

print "Done.\n";
  
