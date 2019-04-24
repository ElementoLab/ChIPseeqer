#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;




my $rscript = "m <- read.csv(\"$ARGV[0]\", header=T, row.names=1, sep=\"\\t\", check.names=T)\n";
$rscript   .= "print(cor.test(m[,$ARGV[1]], m[,$ARGV[2]], method=\"pearson\"))\n";


open OUT, "| R --slave ";
print OUT $rscript;
close OUT;
