#!/usr/bin/perl

my $t = shift @ARGV;


my %H = ();
my %P = ();
foreach my $f (@ARGV) {
  open IN, $f;
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    if ($a[4] < $t) {
      $H{$a[0]} ++;
    }
    $P{$a[0]} = 1;
  }
  close IN;

}

my $cnt_t = @ARGV / 2;
#if (@ARGV == 1) {
#  $cnt_t = 0;
#}
foreach my $p (keys(%P)) {
  if ($H{$p} > $cnt_t) {
    print "$p\t1\n";
  } else {
    print "$p\t0\n";
  }
}
