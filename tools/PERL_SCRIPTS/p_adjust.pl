#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";

my $col = undef;

$col = $ARGV[1] if ($ARGV[1] ne "");

my @p = ();
my $h = <IN>; 
chomp $h;
print "$h\tpadj\n";

my @MAT = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if (!defined($col)) {
    $col = $#a;
  }
  #print "$a[$col]\n"; 
  push @p, $a[$col];
  push @MAT, $l;

}
close IN;


my $padj = Sets::BH_adjust(\@p);


my $i = 0;
foreach my $l (@MAT) {
  print "$l";
  print sprintf("\t%3.2e", $padj->[$i]);
  print "\n";
  $i++;
}


