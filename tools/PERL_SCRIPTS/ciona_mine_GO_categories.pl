BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Sets;
use GroupEnrichment;
use strict;

my $go = GroupEnrichment->new;

$go->setGroups($ARGV[0]);

if ($ARGV[2] != -1) {
    $go->setGroupDesc($ARGV[1]);
}

$go->setMinGroupSize(5);
$go->setMaxGroupSize(1000);
$go->setBonferroni(1);

my $a_ref_ids = $go->getGroupIds();

foreach my $r (@$a_ref_ids) {

 

  my $ref = $go->getGeneGroup($r);
  
  my $gg = $go->getDesc($r);
  #next if ((scalar(@$ref) < 10) || (scalar(@$ref) > 1000));
  next if ($gg !~ /uscle/);

  print $go->getDesc($r); print "\n";
  
  Sets::writeSet($ref, "toto.txt");

  my $todo = "/home/olly/PROGRAMS/CLUSTERS/enrich -clusterfile toto.txt  -kmerfile /home/olly/DATA/KMERS/6mers.txt -kmersize 6 -fastafile u_1000_c_0_masked.seq";

  system($todo);
  
}
