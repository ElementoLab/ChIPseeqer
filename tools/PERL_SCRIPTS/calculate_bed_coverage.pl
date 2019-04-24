use strict;
open IN, $ARGV[0];
my @cov = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  next if ($a[0] !~ /chr/);
  for (my $i=$a[1]; $i<$a[2]; $i++) {
    $cov[$i] ++;
  }
}
close IN;


my $sum = 0;
for (my $i=0; $i<@cov; $i++) {
  if ($cov[$i] > 0) {
    $sum ++;
  }
}

print "$sum\n";
