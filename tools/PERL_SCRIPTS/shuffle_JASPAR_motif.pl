BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use strict;




my $perm = undef;

open IN, $ARGV[0];

my @MATRIX = ();

my $lcnt   = 0;
my $maxlen = 0;

while (my $r = <IN>) {
  
  #print "$r";

  chomp $r;

  my @a = split /\ +/, $r, -1;
  shift @a if ($a[0] eq "");

  #print join("-", @a) . "\n";


  if ($lcnt == 0) {
    my $m = scalar(@a);
    $m --;
    $perm = Sets::shuffle_array([0..$m]);
  }

  #print join("->", @$perm) . "\n";


  for (my $i=0; $i<@$perm; $i++) {
    my $n = $a[$perm->[$i]];
    push @{ $MATRIX[$lcnt] }, $n;
    my $l = length("$n");
    if ($l > $maxlen) {
      $maxlen = $l;
    }
  }

  $lcnt ++;
}

close IN;

$maxlen ++;
foreach my $r (@MATRIX) {
  foreach my $s (@$r) {
    printf("%$maxlen" . "d", $s);
  }
  print "\n";
}
