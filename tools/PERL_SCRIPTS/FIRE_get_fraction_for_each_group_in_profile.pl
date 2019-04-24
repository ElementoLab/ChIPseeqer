#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;
use Getopt::Long;

my $f = $ARGV[0];
my $skip = undef;
my $perc = 1;
GetOptions("skip=s" => \$skip,
	   "perc=s" => \$perc);

my @a_skip = undef;
if (defined($skip)) {
  @a_skip = split /\,/, $skip;
}

open IN, $f or die "Cannot open $f\n";
my $l = <IN>;
my %H = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if (defined($skip) && Sets::in_array($a[1], @a_skip)) {
    next;
  }
  $H{$a[1]} ++; 
}
close IN;


my $tot = 0;
foreach my $k (keys(%H)) {
  $tot += $H{$k};
}

foreach my $k (sort(keys(%H))) {
  my $txt = "%%";
  if ($perc == 0) {
    $txt = "";
  }
  printf("$k\t%3.2f$txt\n", 100*$H{$k}/$tot);
}



