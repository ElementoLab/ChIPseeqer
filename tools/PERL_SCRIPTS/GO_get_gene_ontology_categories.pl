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
$go->setMaxGroupSize((defined($ARGV[5])?$ARGV[5]:500));




my $a_ref_genes = Sets::readSet($ARGV[0]);
foreach my $r (@$a_ref_genes) {

  print "$r\n";

  my $a_ref = $go->getGeneAnnotation($r);
    
  foreach my $s (@$a_ref) {
    print "\t" . $go->getDesc($s);
    print "\n";
  }

  

}


