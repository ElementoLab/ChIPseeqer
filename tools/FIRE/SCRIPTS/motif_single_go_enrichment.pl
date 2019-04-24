use lib "$ENV{FIREDIR}/SCRIPTS";

#
# 1. get output of micombine as input
# 2. do a GO analysis on them
# 3
#

use Sets;
use GroupEnrichment;
use Getopt::Long;
use Table;

use strict;


my $re             = undef;
my $fastafile      = undef;
my $goindex        = undef;
my $gonames        = undef;
my $maxcatsize     = 500;
my $pmin           = 0.01;
my $outgo          = undef;
my $outfile        = undef;

my $profiles       = undef;
my $expfile        = undef;
my $summaryfile    = undef;
my $clusterfile    = undef;
my $usemodule      = undef;

if (@ARGV == 0) {
  die "Usage: perl motif_single_go_enrichment.pl --summaryfile=FILE --expfile=FILE --profiles=FILE --goindex=FILE --gonames=FILE\n";
}

GetOptions ('clusterfile=s'    => \$clusterfile,
	    'summaryfile=s'    => \$summaryfile,
	    'profiles=s'       => \$profiles,
	    'expfile=s'        => \$expfile,
	    'goindex=s'        => \$goindex,
	    'maxcatsize=s'     => \$maxcatsize,
	    'outgo=s'          => \$outgo,
            'outfile=s'        => \$outfile,
	    'pmin=s'           => \$pmin,
	    'usemodule=s'      => \$usemodule,
	    'gonames=s'        => \$gonames);


use Table;

my $ta = Table->new;


#
# read in clusters associated to each motif
#
my %CLU = ();
my @MOTIFS = ();
$ta->loadFile($summaryfile);
my $a_ref_sum = $ta->getArray();
foreach my $r (@$a_ref_sum) {
  push @MOTIFS, $r->[0];
  for (my $i=12; $i<@$r; $i++) {
    $CLU{$r->[0]}{$r->[$i]} = 1;
  }
}


if (defined($clusterfile)) {
  #
  # read in clusters of motifs
  #
  @MOTIFS = ();
  $ta->loadFile($clusterfile);
  my $a_ref_mcl = $ta->getArray();
  foreach my $r (@$a_ref_mcl) {
    push @MOTIFS, $r->[0];
  }

}

#
#  read in expfile
#
$ta->loadFile($expfile);
my $h_ref_exp = $ta->getIndexKV(0,1);
my $a_ref_exp = $ta->getColumn(0);

#
# load profiles
#
my %CLU_GENES = ();

$ta->loadFile($profiles);
my $a_ref_pro = $ta->getArray();
foreach my $r (@$a_ref_pro) {
  # get the expression cluster
  my $c = $h_ref_exp->{ $r->[1] };
  my $m = $r->[0];
  #print "$m\t$c ($r->[1])\n";
  next if (!defined($c));
  if (defined($CLU{$m}{$c})) {
    push @{ $CLU_GENES{ $m }}, $r->[1];
  }
}

if (defined($outfile)) {
  open OUT, ">$outfile" or die "Cannot open $outfile.\n";
}

open OUT2, ">$outgo" or die "cannot open outgo $outgo.\n";


foreach my $k (@MOTIFS) {



  my $s  = Sets::removeDuplicates($CLU_GENES{$k});
  my $nn = scalar(@$s);

  
  my $go = GroupEnrichment->new;
  
  if (defined($usemodule)) {
    $go->setUseModule($usemodule);
  }

  $go->setGeneUniverse($a_ref_exp);
  $go->setGroups($goindex);
  $go->setGroupDesc($gonames);
  $go->setMinGroupSize(5);
  $go->setMaxGroupSize($maxcatsize);
  $go->setBonferroni(1);
    
  my $a_ref = $go->getGroupEnrichment($s, -1, (defined($pmin)?$pmin:0.01));
  
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
    print OUT2 "$k\t$r->[6], $pv\t$r->[5]\n"; 
    print "$k\t$r->[6], $pv\n";  
    
  } else {
    print OUT2 "$k\t\n";
    print "$k\t\n";
  }
  
  if (defined($outfile)) {
    foreach my $r (@$a_ref) {
      my $pv = sprintf("%1.2e", $r->[0]);
      print OUT "$k\t$pv\t$r->[1]\t$r->[2]\t$r->[3]\t$r->[4]\t$r->[6]\n";
    }
  }
  
}


close OUT if (defined($outfile));
close OUT2 if (defined($outgo));




