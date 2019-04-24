#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;
use strict;

if (@ARGV == 0) {
  die "args: file numcol1 numcol2\n";
}
my $fdr = 0.1;
my $fold = 2;
my $l2fold = Sets::log2($fold);
my $fi  = $ARGV[0];

my $x1  = 1;
my $x2  = $ARGV[1];
my $y1  = $x2+1;
my $y2  = $ARGV[1] + $ARGV[2];

my $tt = "two.sided";
if ($ARGV[3]) {
  $tt = $ARGV[3];
}

my $fo = "$fi.t.test.$ARGV[1].$tt.$ARGV[2].fdr$fdr";

my $txt = "
m <- read.csv(\"$fi\", sep=\"\\t\", row.names=1, header=T)
pv <- as.matrix(t(apply(m, 1, function(x) { tt <- t.test(x[$x1:$x2], x[$y1:$y2], alternative=\"$tt\" ); c(tt\$p.value, tt\$statistic) } )))
logfold <- log2( rowMeans(m[,$y1:$y2]) / rowMeans(m[,$x1:$x2]) )

idx   <- rep(-1,dim(m)[1])
padj  <- p.adjust(pv[,1], method=\"BH\")
idxdo <- padj<$fdr & pv[,2]>0
idxup <- padj<$fdr & pv[,2]<0
idxst <- padj>$fdr

idxdo <- idxdo & (logfold < -$l2fold)
idxup <- idxup & (logfold > $l2fold)

idx[idxdo] <- 2  # decrease
idx[idxup] <- 1  # increase
idx[idxst] <- 0

idx <- as.matrix(idx)
rownames(idx) <- rownames(m)
write.table(idx, file=\"$fo\", sep=\"\\t\", row.names=T, col.names=NA, quote=F)

idx <- rep(-1,dim(m)[1])
idx[idxdo] <- 1  # decrease
idx[idxup] <- 0  # increase
idx[idxst] <- 0

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


