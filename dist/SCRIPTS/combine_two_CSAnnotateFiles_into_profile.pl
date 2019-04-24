#!/usr/bin/perl


use lib "$ENV{CHIPSEEQERDIR}";

use Table;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args --genefile1=FILE --genefile2=FILE --geneparts=STR\n";
}
my $genefile1 = undef;
my $genefile2 = undef;
my $merge     = "AND";
my $geneparts = "P";

GetOptions("genefile1=s" => \$genefile1,
           "genefile2=s" => \$genefile2,
	   "geneparts=s" => \$geneparts);

my $ta = Table->new;

# load first file
$ta->loadFile($genefile1);
my $h = $ta->shift;
my $idxgp = -1;
my $j = 0;
foreach my $r (@$h) {
  if ($r eq $geneparts) {
    $idxgp = $j;
  }
  $j++;
}
if ($idxgp == -1) {
  die "Could not find gene part $geneparts\n";
}

my $a_ref = $ta->getArray();

# load second file
$ta->loadFile($genefile2);
my $h_ref = $ta->getIndexKV(0,$idxgp);

shift @$a_ref;
print "GENE\tEXP\n";
foreach my $r (@$a_ref) {

  my $i = 0;
  if (($r->[$idxgp] == 0) && ($h_ref->{ $r->[0] } == 0)) {
    $i = 0;
  } elsif (($r->[$idxgp] == 0) && ($h_ref->{ $r->[0] } == 1)) {
    $i = 1;
  } elsif (($r->[$idxgp] == 1) && ($h_ref->{ $r->[0] } == 0)) {
    $i = 2;
  } elsif (($r->[$idxgp] == 1) && ($h_ref->{ $r->[0] } == 1)) {
    $i = 3;
  }
  print "$r->[0]\t$i\n";

}

