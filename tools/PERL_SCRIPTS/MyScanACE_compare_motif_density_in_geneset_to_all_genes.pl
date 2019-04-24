#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use strict;
use Statistics::Test::WilcoxonRankSum;




if (@ARGV == 0) {
  die "Args: scanaceprofile geneset\n";
}

# load geneset
my %GENES = ();
open IN, $ARGV[1] or die "Cannot open $ARGV[1]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $GENES{$a[0]} = 1;
}
close IN;


# go thru profile
my @gs_dens = ();
my @ag_dens = ();

open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
my $l = <IN>;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if (defined($GENES{$a[0]})) {
    push @gs_dens, $a[1];
  } else {
    push @ag_dens, $a[1];
  }
}
close IN;


my $wilcox_test = Statistics::Test::WilcoxonRankSum->new;

my $gs_avg = Sets::average(\@gs_dens);
my $ag_avg = Sets::average(\@ag_dens);

my $pf = undef;
my $st = undef;
if (($gs_avg == 0) || ($ag_avg == 0)) {
  #print "$ARGV[0]\tNA\n";
  $pf = "NA";
  #exit;
} else {
  
  #print join(" ", @ag_dens) . "\n";
  
  $wilcox_test->load_data(\@gs_dens, \@ag_dens);
  
  my $prob = $wilcox_test->probability();

  $pf = sprintf '%5.3e', $prob; # prints 0.091022

  if (($gs_avg > $ag_avg) && ($prob < 0.01/25)) {
    $st = "***";
  }
}

$gs_avg = sprintf("%3.2f", $gs_avg);
$ag_avg = sprintf("%3.2f", $ag_avg);

my $fs = Sets::filename($ARGV[0]);
$fs =~ s/\_pwm.*$//;

print "$fs\t$gs_avg\t$ag_avg\t$pf\t $st\n";

#print $wilcox_test->probability_status();
