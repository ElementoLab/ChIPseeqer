my $f1 = $ARGV[0];
my $f2 = $ARGV[1];

my $s1 = $f1; $s1 =~ s/\_master.*$//;
my $s2 = $f2; $s2 =~ s/\_master.*$//;

my $cm = "count_cmp_$s1\_$s2.txt";

my $todo = undef;
$todo = "expression_concatenate_matrices.pl $f1 $f2 > $cm";
system($todo) == 0 or die "Cannot exec $todo\n";


my $Rscript = "
m <- read.csv(\"$cm\", sep=\"\t\", header=T, row.names=1)
e <- 0.0001;
foldChange <- (m[,3]+e)/(m[,6]+e)
cmpFreq <- function(y, n) { x <- as.numeric(y); if (is.na(x[1]) | is.na(x[4])) { 1 } else { p <- prop.test(c(x[1], x[4]), c(x[2], x[5]), alternative=\"greater\")\$p.value; min(1,(p*n)) } }
pv <- apply(m, 1, cmpFreq, dim(m)[1])
mb <- cbind(m, foldChange, pv)
write.table(mb[order(mb[,8]),], file=\"$cm.stats.txt\", sep=\"\t\", quote=F, row.names=T, col.names=NA)
";

open OUT, "| R --slave " or die "Cannot open R this way\n";
print OUT $Rscript;
close OUT;

# cleanup outfile

open IN, "$cm.stats.txt";
my $l = <IN>; chomp $l;
my @a = split /\t/, $l, -1;
$a[0] = "RepeatFamily";
$l = join("\t", @a);

open OUT, ">$cm.stats.formatted.txt";

print "$l\tSignif\n";
print OUT "$l\tSignif\n";

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $a[7] = sprintf("%3.2f", $a[7]);
  $a[8] = sprintf("%3.2e", $a[8]);
  if ($a[8] <= 0.05) {
    push @a, "***";
  } else {
    push @a, "   ";
  }
  print join("\t", @a) . "\n";
  print OUT join("\t", @a) . "\n";
}

close OUT;
close IN;
