#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

use Getopt::Long;

#GetOptions("matrix=s"   => $matrixfile,

if (@ARGV == 0) {
  die "Args: matrix phenotypes\n";
}

#
# load matrix
#
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
$ta->processHeader();




#
# load phenotypes
#
$ta->setPhenotypeTable($ARGV[1]);

#
# reaorganize using phenotypes
#
if ($ARGV[2] eq "") {
  $ta->reorderColumnsUsingPhenotype;
} else {
  my @a = split /\,/, @ARGV[2];
  $ta->reorderColumnsUsingPhenotype(\@a);
}
#
# print new matrix
#
my $a_ref = $ta->getArray();
my $a_ref_h = $ta->getHeader();
print join("\t", @$a_ref_h) . "\n";

foreach my $r (@$a_ref) {
  print join("\t", @$r) . "\n";
  
}

