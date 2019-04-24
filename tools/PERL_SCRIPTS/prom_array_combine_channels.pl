use strict;

if (@ARGV == 0) {
  die "Args: suffix (will add _532 and _635)\n";
}

my $suf = $ARGV[0];

#
# load 532 data
#
my $f532 = "$suf\_532.pair";
if (! -e $f532) {
  $f532 .= ".txt";
}

#
# load 635 data
#
my $f635 = "$suf\_635.pair";
if (! -e $f635) {
  $f635 .= ".txt";
}

my %H1 = ();
open IN, $f635;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;
  if ($a[9] == 0.0) {
    $a[9] = 0.1;
  }

  $H1{$a[3]} = $a[9];
}
close IN;

open IN, $f532;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;

  #print "$a[3]\t" . $H1{$a[3]} . "\t$a[9]\n";

  if ($a[9] == 0.0) {
    $a[9] = 0.1;
  }

  my $rat = sprintf("%4.3f", log($H1{$a[3]} / $a[9]));
  print "$a[3]\t$rat\n";
  
}
close IN;

