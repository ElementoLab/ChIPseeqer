BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use Sets;
use GroupEnrichment;
use strict;
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref_lines = $ta->getArray();


my $go = GroupEnrichment->new;

$go->setGroups($ARGV[1]);

if ($ARGV[2] != -1) {
    $go->setGroupDesc($ARGV[2]);
}

$go->setMinGroupSize(5);
$go->setMaxGroupSize((defined($ARGV[5])?$ARGV[5]:500));
$go->setBonferroni(1);

#$go->setVerbose(1);



die("define N ..\n") if !defined($ARGV[3]);

my $line = 0;
foreach my $a_ref_genes (@$a_ref_lines) {

  my $a_ref_genes_unique = Sets::removeDuplicates($a_ref_genes);

  #print sprintf("%10s%4s%4s  %s\n", "pvalue", "cnt", "exp", "category");    
  my $a_ref = $go->getGroupEnrichment($a_ref_genes_unique, $ARGV[3], (defined($ARGV[4])?$ARGV[4]:0.00001));
  foreach my $r (@$a_ref) {
    my $cnt = $r->[1];
    my $exp = $r->[3] * $r->[2] / $r->[4];
  
    my $pv  = sprintf("%1.2e", $r->[0]);
    my ($e) = $pv  =~ /e\-(\d{2})/;
    
    
    $e = int($e)-1;
    
    
    #$pv = eval("10**$e");
    
    #$pv = 1.0/$pv;
    
    #if (int($e) > 3) {
    #  $pv = sprintf("<%2.2e", $pv);
    #} else {
    #  $pv = sprintf("<%s", $pv);
    #}
    
    print sprintf("%s\t%d\t%3.2e\t%d\t%d\t%d\t%s\n", $ARGV[0], $line, $pv, $r->[2], $cnt, $exp, $r->[6]);    
    #print sprintf("pv=%3.2e\tov=%d\tca=%d\tse=%d\tN=%d\ttxt=%s\t%s\n", $r->[0], $r->[1], $r->[2], $r->[3], $r->[4], $r->[5], $r->[6]);
  }

  $line ++;
}
  
