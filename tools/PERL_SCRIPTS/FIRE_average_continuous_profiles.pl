use strict;

if (@ARGV == 0) {
  die "Args: f1 f2\n";
}

my $f1 = $ARGV[0];
my $f2 = $ARGV[1];

my %H1 = ();
open IN, $f1;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;
  $H1{$a[0]} = $a[1];
}
close IN;

open IN, $f2;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;

  my $rat = ($H1{$a[0]} + $a[1]) / 2; 
  print "$a[0]\t$rat\n";
  
}
close IN;

