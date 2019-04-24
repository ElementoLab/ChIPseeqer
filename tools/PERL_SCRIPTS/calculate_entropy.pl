use lib qw(/home/elemento/PERL_MODULES);
use Sets;

while (my $l = <STDIN>) {
  chomp $l;
  my @a = split /\t/, $l;

  my $n = shift @a;

  my $s = Sets::stddev(\@a);

  print "$n\t$s\n";
}
