BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;

open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";

my $l = <IN>;
print "GENE\tVAR\n";

while (my $l = <IN>) {
  chomp $l;

  my @b = split /\t/, $l, -1;
  my $n = shift @b;

  my $a = Sets::average(\@b);
  my $s = Sets::stddev (\@b);

  if (abs($a) > 1e-4) {
    my $o = sprintf("%4.3f", $s / $a);
    print "$n\t$o\n";
  }

}

close IN;
