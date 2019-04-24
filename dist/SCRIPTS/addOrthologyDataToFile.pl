#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args --file=FILE --orthology=FILE --do=STR[add|replace]\n";
}
my $orthology = "$ENV{CHIPSEEQERDIR}/DATA/hg18/hum_mou_NM_orthologs_29NOV20009_20APR2010.txt";
my $file      = undef;
my $do        = "add";
GetOptions("file=s"      => \$file,
	   "do=s"        => \$do,
           "orthology=s" => \$orthology);

# load orthol
open IN, $orthology or die "Cannot open $orthology\n";
#print STDERR "# using $orthogy\n";
my %H = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if ($a[1] !~ /not/) {
    $H{$a[1]} = $a[0];
  }
}
close IN;

# file
open IN, $file or die "Cannot open $file\n";
my $l = <IN>; chomp $l;
if ($do eq "add") {
  print "$l\tORTH\n";
} else {
  my @a = split /\t/, $l, -1;
  $a[0] = "ORTH";
  $l = join("\t", @a);
  print "$l\n";
}

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  #print $l;
  if (defined($H{$a[0]})) {
    if ($do eq "add") {
      print "$l\t$H{$a[0]}\n";
    } else {
      $a[0] = $H{$a[0]};
      $l = join("\t", @a);
      print "$l\n";
    }
  } else {
   # print "H{$a[0]} = $H{$a[0]} not found?\n";
  } 

}
close IN;


