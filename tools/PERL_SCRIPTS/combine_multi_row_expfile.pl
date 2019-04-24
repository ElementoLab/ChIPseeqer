use strict;

open IN, $ARGV[0] or die "Cannot open file";
my $l = <IN>;
my $max1 = -1;
my $max2 = -1;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if ($a[1] > $max1) {
    $max1 = $a[1];
  }
  if ($a[2] > $max2) {
    $max2 = $a[2];
  }
}
close IN;

my %H = ();
my $idx = 0;
for (my $i=0; $i<=$max1; $i++) {
  for (my $j=0; $j<=$max2; $j++) {
    print STDERR "$i-$j\t$idx\n";
    $H{"$i-$j"} = $idx;
    $idx ++;
  }
}


open IN, $ARGV[0] or die "Cannot open file";
my $l = <IN>;
chomp $l;
my @a = split /\t/, $l, -1;
print "$a[0]\tEXP\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $k = "$a[1]-$a[2]";
  print "$a[0]\t$H{$k}\n";
}
close IN;

