BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

open IN, $ARGV[0];
my $i = 0;
while (my $l = <IN>) {


  if ($i == $ARGV[1]) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    print "$a[$ARGV[2]]\n";
  }
  
  $i++;
}
close IN;
