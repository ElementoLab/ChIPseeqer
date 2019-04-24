#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


if (@ARGV == 0) {
  die "Args: dir lane\n";
}

my $dir  = $ARGV[0];
my $lane = $ARGV[1];

my $f1 = "$dir/s_$lane\_1_sequence.txt.gz";
my $maxlines = 100000;

my %B = ();
if (-e $f1) {
  my $cnt = 0;
  open IN, "zcat $f1 | ";
  while (my $l = <IN>) {
    next if ($l !~ /\@/);
    chomp $l;
    my @a = split /\t/, $l, -1;
    my ($bc) = $a[0] =~ /\#(.+?)\//;
    next if ($bc =~ /N/);
    $B{$bc}++;
    $cnt++;
    last if ($cnt == $maxlines);
  }
  close IN;

  my @k = keys(%B);

  foreach my $c (@k) {
    print "$c\t$B{$c}\n";
  }
  
} else {
  die "Could not find $f1\n";
}


