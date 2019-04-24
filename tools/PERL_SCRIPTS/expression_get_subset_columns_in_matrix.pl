BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

my $a_ref_cols = Sets::readSet($ARGV[0]);

open IN, $ARGV[1] or die "Cannot open $ARGV[1]\n";
while (my $l = <IN>) {

  my @a = split /\t/, $l, -1;
  my $n = shift @a;

  my @b = ();
  push @b, $n;
  foreach my $i (@$a_ref_cols) {
    push @b, $a[$i];
  }
  
  print join("\t", @b) . "\n";

}
close IN;
