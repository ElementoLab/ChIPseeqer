#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


my $file = shift @ARGV;
open IN, $file;
my $l = <IN>;
chomp $l;
my @a = split /\t/, $l;
my $r1 = shift @a;


my $i = 0;

my @COL_SETS = ();
foreach my $r (@a) {
  my $j = 0;
  foreach my $re (@ARGV) {
    if ($r =~ /$re/)  {
      push @{ $COL_SETS[ $j ] }, $i;      
    }
    $j++;
  }
  $i++;
}

print "$r1";
foreach my $cs (@COL_SETS) {
  foreach my $i (@$cs) {
    print "\t" . $a[$i];
  }
}
print "\n";


while ($l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $r1 = shift @a;
  print "$r1";
  foreach my $cs (@COL_SETS) {
    foreach my $i (@$cs) {
      print "\t" . $a[$i];
    }
  }
  print "\n";
}





