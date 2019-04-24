BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use strict;




my $perm = undef;

open IN, $ARGV[0];

my @MATRIX = ();

my $lcnt   = 0;
my $maxlen = 0;

my $l = <IN>;
print $l;

while (my $r = <IN>) {
  
  #print "$r";

  chomp $r;

  my @a = split /\t/, $r, -1;
  
  my $nt = shift @a;

  if ($lcnt == 0) {
    my $m = scalar(@a);
    $m --;
    $perm = Sets::shuffle_array([0..$m]);
  }

  print "$nt";
  for (my $i=0; $i<@$perm; $i++) {
    my $n = $a[$perm->[$i]];
    print "\t$n";
  }
  print "\n";
  $lcnt ++;
}

close IN;
