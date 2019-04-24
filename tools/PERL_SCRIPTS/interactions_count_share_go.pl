BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use Sets;
use GroupEnrichment;
use strict;

my $go = GroupEnrichment->new;

$go->setGroups($ARGV[1]);

if ($ARGV[2] != -1) {
    $go->setGroupDesc($ARGV[2]);
}

#$go->setMinGroupSize(5);
$go->setMaxGroupSize(250);

open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;

  #print "$a[0]\t$a[1]\n";

  my $a_ref1 = $go->getGeneAnnotation($a[0]);
  my $a_ref2 = $go->getGeneAnnotation($a[1]);
  
  next if (!$a_ref1 || !$a_ref2);
  
  #print join("\t", @$a_ref1); print "\n";

  my $a_ref3 = Sets::getOverlapSet($a_ref1, $a_ref2);
  
  next if (scalar(@$a_ref3) == 0);

  print "$a[0]\t$a[1]\t"; print join("\t", @$a_ref3); print "\n";

}
close IN;


