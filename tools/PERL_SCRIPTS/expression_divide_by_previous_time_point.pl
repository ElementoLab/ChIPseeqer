#BEGIN{ $home = `echo \$HOME`; chomp $home}
#use lib "$home/PERL_MODULES";
use lib qw(/home/elemento/PERL_MODULES);

use Sets;
use strict;

open IN, $ARGV[0];
my $l = <IN>;
chomp $l;
my @a = split /\t/, $l, -1;
my $c = shift @a;
shift @a;
print "$c\t" . join("\t", @a) . "\n";

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $c = shift @a;
  
  my @b = ();
  my $n = scalar(@a);
  for (my $i=1; $i<$n; $i++) {
    my $o = sprintf("%4.3f", log( $a[ $i ] / $a[ $i - 1 ] ) / log(2.0));
    push @b, $o;
  }
  
  print "$c\t" . join("\t", @b) . "\n";
}
 
close IN;

   
