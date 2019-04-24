BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

open IN, $ARGV[0];
my $l = <IN>;
print "$l";



while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $r = \@a;
  my $n = shift @$r;
  
  my $max = -100000;
  foreach my $s (@$r) {
    if ($s > $max) {
      $max = $s;
    }    
  }
  
  if ($max >= $ARGV[1]) {
    print "$n\t" . join("\t", @$r) . "\n";
  }

}
close IN;


