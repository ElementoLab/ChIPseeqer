#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

if (@ARGV == 0) {
  die "Args : dataset gene\n";
}	

my $txt = "
m <- read.csv(\"$ARGV[0]\", sep=\"\\t\", row.names=1, header=T, check.names=F)
v <- apply(m, 1, function(x) { cc<-cor.test(as.numeric(m[\"$ARGV[1]\",]),as.numeric(x)); c(cc\$estimate, cc\$p.value) })
v <- t(v)
v <- cbind(v, p.adjust(as.numeric(v[,2]), method=\"BH\"))
ve <- v[,3]<0.1 & v[,1]>0
ve[ve == F] <- 0
write.table(ve, file=\"$ARGV[0].Rcorrel_with_$ARGV[1].fdr0.1.pos\", sep=\"\t\", quote=F, row.names=T, col.names=NA)

"; 

my $tmpfile = Sets::getTempFile("/tmp/Rscript");
Sets::writeText($txt, $tmpfile);

system("R CMD BATCH $tmpfile") == 0 or die "Cannot exec R script ?\n";

#system("cat $file.txt");

unlink $tmpfile;
