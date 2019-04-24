#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args: matrix phenotypes label1 label2 \n";
}

my $ta = Table->new;
$ta->loadFile($ARGV[0], 0, 1);
#$ta->getHeader();
$ta->setPhenotypeTable($ARGV[1]);
my $a_ref = $ta->getArray();

print "GENE\tT\n";

foreach my $r (@$a_ref) {

  #my $a_ref_c1     = $ta->getValuesForPhenotype($r, "GC B-Like");
  #my $a_ref_c2     = $ta->getValuesForPhenotype($r, "Activated B-like");
 
  my $a_ref_c1     = $ta->getValuesForPhenotype($r, $ARGV[2]);
  my $a_ref_c2     = $ta->getValuesForPhenotype($r, $ARGV[3]);

  my $n1           = @$a_ref_c1;
  my $n2           = @$a_ref_c2;

  # print STDERR "$n1\t$n2\n";
 
  #my $a_ref_c1_nom = Sets::logArray(Sets::removeMissingValuesInArray($a_ref_c1));
  #my $a_ref_c2_nom = Sets::logArray(Sets::removeMissingValuesInArray($a_ref_c2));

  my $a_ref_c1_nom = Sets::removeMissingValuesInArray($a_ref_c1);
  my $a_ref_c2_nom = Sets::removeMissingValuesInArray($a_ref_c2);


  my $t            = undef;

  if (defined($ARGV[4])) {
    my $m1           = Sets::average($a_ref_c1_nom);
    my $m2           = Sets::average($a_ref_c2_nom);
    $t = sprintf("%4.3f", Sets::log2($m1/$m2));
  } else {
    $t = Sets::t_statistic($a_ref_c1_nom, $a_ref_c2_nom);
  }


  print "$r->[0]\t$t\n";
  
  #my $a_ref_c3 = $ta->getPhenotypeColIndices("B-[Ll]ike");
  #print join("\t", @$a_ref_c3) . "\n";

  #print join("\t", @$a_ref_c1) . "\n";
  #print join("\t", @$a_ref_c2) . "\n";
  
  #exit;

}

