#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

my %H = ();
open IN, $ARGV[1];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $H{$a[0]} = 1;

}
close IN;


open IN, $ARGV[0];
my %INT = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $INT{$a[0]} ++;
  $INT{$a[2]} ++;
}
close IN;


foreach my $k (keys(%H)) {
  
  if (defined($INT{$k}) && ($INT{$k} == 1)) {
    print "$k\t1\n";
  } else {
    print "$k\t0\n";
  }	

}
