open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";

my $l = <IN>; chomp $l;
my @a = split /\t/, $l;
shift @a;
print join("\t", @a) . "\n";
my @COLSUMS = ();

my $n = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  shift @a;
  for (my $i=0; $i<@a; $i++) {
    $COLSUMS[$i] += $a[$i];
  }
  $n++;
}
close IN;

foreach my $r (@COLSUMS) {
  $r = sprintf("%3.2f", $r/$n);
}

print join("\t", @COLSUMS) . "\n";

