BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use Sets;
use GroupEnrichment;
use Table;
use strict;


my $go = GroupEnrichment->new;

$go->setGroups($ARGV[1]);

if ($ARGV[2] != -1) {
    $go->setGroupDesc($ARGV[2]);
}

$go->setMinGroupSize(5);
$go->setMaxGroupSize(1000);
#$go->setBonferroni(1);

#$go->setVerbose(1);


die("define N ..\n") if !defined($ARGV[3]);

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my @a_sets = ();
foreach my $r (@$a_ref) {
  push @{ $a_sets[ $r->[1] ] }, $r->[0];
}

#
#  
#
my $cnt = 0;
foreach my $s (@a_sets) {
  print "Cluster $cnt, "; print scalar(@$s); print " members\n"; $cnt++;
  my $a_ref = $go->getGroupEnrichment($s, $ARGV[3], 0.01);
  foreach my $r (@$a_ref) {
    print sprintf("pv=%3.2e\tov=%d\tca=%d\tse=%d\tN=%d\ttxt=%s\t%s\n", $r->[0], $r->[1], $r->[2], $r->[3], $r->[4], $r->[5], $r->[6]);
  }
  #<STDIN>;
}


