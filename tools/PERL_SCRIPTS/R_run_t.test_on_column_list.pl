#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;
use strict;

if (@ARGV == 0) {
  die "args: file collist [ tail ]\n";
}
my $fdr = 0.01;
my $fi  = $ARGV[0];

my $fc  = $ARGV[1];
open IN, $fc or die "Cannot open $fc\n";
my %CO = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @{$CO{$a[1]}}, "\"$a[0]\"";
}
close IN;


my $c0 = "c(" . join(",", @{$CO{"0"}}) . ")";
my $c1 = "c(" . join(",", @{$CO{"1"}}) . ")";

my $tt = "two.sided";
if ($ARGV[2]) {
  $tt = $ARGV[2];
}

my $ffc = Sets::filename($fc);

my $fo = "$fi.t.test.$ffc.fdr$fdr";

my $txt = "
m <- read.csv(\"$fi\", sep=\"\\t\", row.names=1, header=T)
pv <- as.matrix(t(apply(m, 1, function(x) { tt <- t.test(x[$c0], x[$c1], alternative=\"$tt\" ); c(tt\$p.value, tt\$statistic) } )))

idx <- rep(-1,dim(m)[1])
idx[p.adjust(pv[,1], method=\"BH\")<$fdr & pv[,2]>0] <- 2  # decrease
idx[p.adjust(pv[,1], method=\"BH\")<$fdr & pv[,2]<0] <- 1  # increase
idx[p.adjust(pv[,1], method=\"BH\")>$fdr] <- 0
idx <- as.matrix(idx)
rownames(idx) <- rownames(m)
write.table(idx, file=\"$fo\", sep=\"\\t\", row.names=T, col.names=NA, quote=F)

idx <- rep(-1,dim(m)[1])
idx[p.adjust(pv[,1], method=\"BH\")<$fdr & pv[,2]>0] <- 1  # decrease
idx[p.adjust(pv[,1], method=\"BH\")<$fdr & pv[,2]<0] <- 0  # increase
idx[p.adjust(pv[,1], method=\"BH\")>$fdr] <- 0
idx <- as.matrix(idx)
rownames(idx) <- rownames(m)
write.table(idx, file=\"$fo.down.txt\", sep=\"\\t\", row.names=T, col.names=NA, quote=F)

idx <- rep(-1,dim(m)[1])
idx[p.adjust(pv[,1], method=\"BH\")<$fdr & pv[,2]>0] <- 0  # decrease
idx[p.adjust(pv[,1], method=\"BH\")<$fdr & pv[,2]<0] <- 1  # increase
idx[p.adjust(pv[,1], method=\"BH\")>$fdr] <- 0
idx <- as.matrix(idx)
rownames(idx) <- rownames(m)
write.table(idx, file=\"$fo.up.txt\", sep=\"\\t\", row.names=T, col.names=NA, quote=F)




";

my $ft = Sets::getTempFile("/tmp/Rscript");
Sets::writeText($txt, $ft);

system("R CMD BATCH $ft") == 0 or die "Cannot exec R script ?\n";

system("cat $fo");

#unlink $ft;


