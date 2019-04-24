use lib "$ENV{FIREDIR}/SCRIPTS";

#
# 1. get genes that have a given motif
# 2. do a GO analysis on them
#
use Sets;
use GroupEnrichment;
use Table;
use Getopt::Long;



use strict;

if (@ARGV == 0) {
  die "USage: perl clusters_go_enrichment.pl --clusters --goindex --gonames --N --maxcatsize --pmin --outfile\n";
}

my $clusters= undef;
my $index   = undef;
my $names   = undef;
my $N       = -1;
my $maxcatsize = 500;
my $pmin    = 0.01;
my $outfile = undef;
my $outgo   = undef;
my $usemodule = undef;

GetOptions ('clusters=s'        => \$clusters,
	    'goindex=s'         => \$index,
	    'gonames=s'         => \$names,
	    'N=s'               => \$N,
	    'maxcatsize=s'      => \$maxcatsize,
	    'pmin=s'            => \$pmin,
	    'outgo=s'           => \$outgo,
	    'usemodule=s'       => \$usemodule,
            'outfile=s'         => \$outfile);



my $ta = Table->new;
$ta->loadFile($clusters);
my $a_ref = $ta->getArray();
shift @$a_ref;

my %GROUPS        = ();
my @gene_universe = ();

foreach my $r (@$a_ref) {
  push @{ $GROUPS{ $r->[1] } }, $r->[0];
  push @gene_universe, $r->[0];
}

if (defined($outfile)) {
  open OUT, ">$outfile" or die "Cannot open $outfile.\n";
}

open OUT2, ">$outgo" or die "cannot open outgo $outgo.\n";

#
#  count
#

#  Sets::printSet($a_ref_genes);

foreach my $gk (sort {$a <=> $b} (keys(%GROUPS))) {

  my $a_ref_genes = $GROUPS{ $gk };

  my $go = GroupEnrichment->new;

  if (defined($usemodule)) {
    $go->setUseModule($usemodule);
  }
    
  $go->setGeneUniverse(\@gene_universe);
  $go->setGroups($index);
  $go->setGroupDesc($names);
  $go->setMinGroupSize(5);
  $go->setMaxGroupSize((defined($maxcatsize)?$maxcatsize:500));
  $go->setBonferroni(1);
  
  die("define N ..\n") if !defined($N);
  
  
  my $a_ref = $go->getGroupEnrichment($a_ref_genes, $N, (defined($pmin)?$pmin:0.00001));

  my $i = 0; 
  if (defined($a_ref->[0]) && ($a_ref->[0]->[2] > $maxcatsize)) {
    # try to find a category of smaller size
    for (my $j=1; $j<@$a_ref; $j++) {
      if ($a_ref->[$j]->[2] <= 500) {
	$i = $j; last;
      }
    }
  }
  
  if (my $r = $a_ref->[$i]) {
    
    my $cnt = $r->[1];
    my $exp = $r->[3] * $r->[2] / $r->[4];
    
    my $pv  = sprintf("%1.2e", $r->[0]);
    my ($e) = $pv  =~ /e\-(\d{2})/;
  
  
    $e = int($e)-1;
    
    
    $pv = eval("10**$e");
  
    $pv = 1.0/$pv;
    
    if (int($e) > 3) {
      $pv = sprintf("p<%1.0e", $pv);
    } else {
      $pv = sprintf("p<%s", $pv);
    }
    
    $r->[6] =~ s/\ +$//;
    print OUT2 "$gk\t$r->[6], $pv\n"; 
    print "$gk\t$r->[6], $pv\n";  
  
  } else {
    print OUT2 "$gk\t\n";
    print "$gk\t\n";
  }
  
  if (defined($outfile)) {
    foreach my $r (@$a_ref) {
      my $pv = sprintf("%1.2e", $r->[0]);
      print OUT "$gk\t$pv\t$r->[1]\t$r->[2]\t$r->[3]\t$r->[4]\t$r->[6]\n";
    }
  }

}

print "$outgo generated.\n";

close OUT if (defined($outfile));
close OUT2 if (defined($outgo));
