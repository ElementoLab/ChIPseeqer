#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

open IN, $ARGV[0];
my @avg = ();

my $l = <IN>;
print $l;

my $i = 0;
my $cntg = 0;
my $b = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  my $n = shift @a;

  for (my $i=0; $i<@a; $i++) {
    $avg[$i] += $a[$i];
  }
  $cntg ++;

  if ((($i % 100) == 0) && ($i != 0)) {
    
    print "$b"; $b++;
    for (my $i=0; $i<@a; $i++) {
      printf("\t%4.3f", $avg[$i] / $cntg);
    }
    print "\n";
    $cntg = 0;
    @avg = ();
    
  }
  
  $i ++;
}
close IN;

