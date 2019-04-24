#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";


use Sets;
use Getopt::Long;
use strict;

my $matches  = undef;
my $genelist = undef;
my $exclude  = undef;

GetOptions("matches=s"  => \$matches,
	   "exclude=s"  => \$exclude,
           "genelist=s" => \$genelist);



my %EXCL = ();
if (defined($exclude)) {
  open IN, $exclude or die "no $exclude\n";
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    $EXCL{$a[0]} = 1;
  }
  close IN;
}

# 
my %GENES = ();
open IN, $matches or die "Cannot open $matches\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $GENES{$a[0]} ++;
}
close IN;

my @v = ();


open IN, $genelist or die "Cannot open $genelist\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if (!defined($GENES{$a[0]})) {
    $GENES{$a[0]} = 0;
  }
  if (defined($exclude) && defined($EXCL{$a[0]})) {
    next;
  }
  push @v, $GENES{$a[0]};
}
close IN;


my $f = Sets::filename($matches);
$f =~ s/\_pwm.+$//;
my $a = sprintf("%3.1f", Sets::average(\@v));
print "$f\t$a\n";

