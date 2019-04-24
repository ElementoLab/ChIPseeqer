BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

open IN, $ARGV[0];
my $l = <IN>; print $l;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  foreach my $r (@a) {
    if ($r eq "") {
      $r = "nan";
    }
  }
  
  print join("\t", @a); print "\n";
  
}
close IN;
