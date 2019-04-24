if (@ARGV == 0) {
  die "Args file col1=control col2=exp (col2-col1)\n";
}

my $col1 = $ARGV[1];
my $col2 = $ARGV[2];

open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
my $l = <IN>;
my @a = split /\t/, $l, -1;

print "$a[0]\t$a[$col2]-$a[$col1]\n";

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $d = $a[$col2]-$a[$col1];
  print "$a[0]\t$d\n";
}
close IN;

