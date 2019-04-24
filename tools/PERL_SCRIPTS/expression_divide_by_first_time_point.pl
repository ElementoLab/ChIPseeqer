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
  my $m = shift @a;
  foreach my $r (@a) {
    $r = sprintf("%4.3f", log($r/$m)/log(2));
  }
  print "$c\t" . join("\t", @a) . "\n";
}
 
close IN;

   
