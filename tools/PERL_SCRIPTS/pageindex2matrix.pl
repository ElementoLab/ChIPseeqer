#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

my %GENES = ();
my %SETS  = ();

open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  my $g = shift @a;

  foreach my $r (@a) {
    $SETS{$r} = 1;
    $GENES{$g}{$r} = 1;
  }


}
close IN;


my @allsets = keys(%SETS);

foreach my $s (@allsets) {
  print "\t$s";
}
print "\n";
foreach my $g (keys(%GENES)) {
  print "$g";
  foreach my $s (@allsets) {
    if ($GENES{$g}{$s} == 1) {
      print "\t1";
    } else {
      print "\t0";
    }
  }
  print "\n";
}
