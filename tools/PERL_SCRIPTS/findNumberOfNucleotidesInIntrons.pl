#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

my $totalchr = 0;
my %CHRLEN = ();
open IN, $ARGV[1];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $totalchr += $a[1];
  $CHRLEN{$a[0]} = $a[1];
}
close IN;



my %CHR = ();
open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if ($a[0] =~ /^I/) {
    push @{$CHR{$a[1]}}, [$a[2],$a[3]];
  }
}
close IN;

my $totalcount = 0;
foreach my $c (keys(%CHR)) {
  
  print STDERR "$c\n";

  next if (!defined($CHRLEN{$c}));
  
  my @seg = @{$CHR{$c}};

  #@seg = sort { $a->[0] <=> $b->[0] } @seg;

  
  my $a_ref_merged = Sets::mergeOverlappingIntervals(\@seg); #Sets::assemble_overlapping_fragments(\@seg);

  my $cnt = 0 ;
  foreach my $r (@$a_ref_merged) {
    $cnt += ($r->[1] - $r->[0]);
  }
  
  my $rat = 100* $cnt/$CHRLEN{$c};
  print "$cnt bp in introns of $c, or $rat%\n";
  $totalcount += $cnt;
  
}

my $trat = 100*$totalcount/$totalchr;
print "$totalcount bp in introns of genome, or $trat% or genome\n";

