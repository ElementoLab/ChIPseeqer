#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;
use strict;

if (@ARGV == 0) {
  die "args: file numcol1 numcol2\n";
}
my $fdr = 0.10;
my $fi  = $ARGV[0];
print STDERR "# fdr = $fdr\n";

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
m <- read.csv(\"$fi\", sep=\"\\t\", row.names=1, header=T, check.names=F)
pv <- as.matrix(t(apply(m, 1, function(x) { if ((sd(x[$x1:$x2]) < 0.0001) | (sd(x[$y1:$y2]) < 0.0001)) { c(1.0, 0.0) } else { tt <- t.test(log2(x[$x1:$x2]+1), log2(x[$y1:$y2]+1), paired=T, alternative=\"$tt\" ); c(tt\$p.value, tt\$statistic) }} )))

fold <- as.matrix(t(apply(m, 1, function(x) { tt <- log2((mean(x[$y1:$y2])+1)/(mean(x[$x1:$x2])+1)); tt } )))
tfold <- as.matrix(t(apply(m, 1, function(x) { tt <- (mean(x[$y1:$y2])+1)/(mean(x[$x1:$x2])+1); tt } )))

idx <- rep(-1,dim(m)[1])
idx[p.adjust(pv[,1], method=\"BH\")<$fdr & pv[,2]>0] <- 2  # decrease
idx[p.adjust(pv[,1], method=\"BH\")<$fdr & pv[,2]<0] <- 1  # increase
idx[p.adjust(pv[,1], method=\"BH\")>$fdr] <- 0
idx <- as.matrix(idx)
rownames(idx) <- rownames(m)
write.table(idx, file=\"$fo\", sep=\"\\t\", row.names=T, col.names=NA, quote=F)

idxtxt <- c();
idxtxt[idx==1] <- \"up\"
idxtxt[idx==2] <- \"down\"
idxtxt[idx==0] <- \"stable\"

padj <- p.adjust(pv[,1], method=\"BH\")
padj[padj > $fdr] <- 0

numcols <- dim(m)[2]
newm <- cbind(m, t=sprintf(\"%3.2f\", pv[,2]), pv=sprintf(\"%4.3f\", pv[,1]), adjpv=sprintf(\"%4.3f\", p.adjust(pv[,1], method=\"BH\")), fold=sprintf(\"%3.2f\", tfold), diffexp.fdr$fdr=idxtxt, idx)

#tmpm <- newm[order(fold, decreasing=T),]
tmpm <- newm[order(pv[,2], decreasing=F),]
newm <- rbind( tmpm[tmpm[,numcols+6]==1,], tmpm[tmpm[,numcols+6]==0,],  tmpm[tmpm[,numcols+6]==2,])[,-(numcols+5)]

write.table(newm, file=\"$fo.fullres.txt\", sep=\"\\t\", row.names=T, col.names=NA, quote=F)


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

system("cat $fo.fullres.txt");

#unlink $ft;


