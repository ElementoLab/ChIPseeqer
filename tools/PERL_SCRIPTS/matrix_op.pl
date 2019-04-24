#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

my $f = shift @ARGV;
my @evals  = ();
my @evalsh = ();
foreach my $e (@ARGV) {
  my ($n1, $op, $n2) = $e =~ /^(\d+)([\/])(\d+)$/, $e;

  my $txt = "Sets::log2((\$a[$n1]+1)$op(\$a[$n2]+1))";
  push @evals, $txt;
  
  $txt = "log(\$a[$n1]$op\$a[$n2])";
  push @evalsh, $txt;

}
open IN, $f or die "Cannot open file $f\n";
my $l = <IN>;
chomp $l;
my @a = split /\t/, $l, -1;
foreach my $e (@evalsh) {
  my $ee = $e;
  $ee =~ s/\//ov/g;		
  eval("print(\"\t$ee\")");
}
print "\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  #print "$a[0]";
  print "$l";
  foreach my $e (@evals) {
    print "\t" . sprintf("%5.4f", eval($e));
  }
  print "\n";

  #print eval("\($a[1]+1)/(\$a[2]+1)");
  
}
close IN;

