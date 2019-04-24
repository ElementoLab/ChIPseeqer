#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


my $mr         = 25;
my $multiplex  = 1;
my $genomesize = 4411532;
my $readlen    = 40;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args --mr=INT --multiplex=INT --genomesize=INT --readlen=INT\n";
}
GetOptions("mr=s"         => \$mr,
           "multiplex=s"  => \$multiplex,
	   "genomesize=s" => \$genomesize,
	   "readlen=s"    => \$readlen);


my $numnt  = $readlen * $mr * 1000000 /$multiplex;

my $cov = int($numnt / $genomesize);


print "$cov" . "X\n";

