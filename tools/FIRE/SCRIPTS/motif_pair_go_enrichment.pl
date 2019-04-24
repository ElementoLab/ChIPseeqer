#
# 1. get output of micombine as input
# 2. do a GO analysis on them
#

use lib "$ENV{FIREDIR}/SCRIPTS";

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
my $pmin           = undef;
my $outgo          = undef;
my $outfile        = undef;
my $micombinefile  = undef;
my $profiles_dna   = undef;
my $profiles_rna   = undef;
my $expfile        = undef;
my $usemodule      = undef;

if (@ARGV == 0) {
  die "Usage: perl motif_pair_go_enrichment.pl --expfile=FILE --pmin=FLT --micombinefile=FILE --profiles=FILE --goindex=FILE --gonames=FILE --maxcatsize=INT\n";
}

GetOptions ('micombinefile=s'  => \$micombinefile,
	    'profiles_dna=s'   => \$profiles_dna,
	    'profiles_rna=s'   => \$profiles_rna,
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
#  read in expfile
#
$ta->loadFile($expfile);
my $a_ref_exp = $ta->getColumn(0);


# load profiles

my %H = ();

if (defined($profiles_dna)) {
  $ta->loadFile($profiles_dna);
  my $a_ref_pro = $ta->getArray();
  foreach my $r (@$a_ref_pro) {
    $H{$r->[0]}{$r->[1]} = 1;
  }
}


if (defined($profiles_rna)) {
  $ta->loadFile($profiles_rna);
  my $a_ref_pro = $ta->getArray();
  foreach my $r (@$a_ref_pro) {
    $H{$r->[0]}{$r->[1]} = 1;
  }
}

#
# load combi file
#
$ta->loadFile($micombinefile);
my $a_ref_combi = $ta->getArray();


if (defined($outfile)) {
  open OUT, ">$outfile" or die "Cannot open $outfile.\n";
}

open OUT2, ">$outgo" or die "cannot open outgo $outgo.\n";


foreach my $s (@$a_ref_combi) {
  if ($s->[3] == 0) {
    my @a1 = keys( %{ $H{ $s->[0] } } );
    my @a2 = keys( %{ $H{ $s->[1] } } );
    my $a3 = Sets::getIntersection(\@a1, \@a2);

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

    
    my $a_ref = $go->getGroupEnrichment($a3, -1, (defined($pmin)?$pmin:0.01));
    
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
      print OUT2 "$s->[0]\t$s->[1]\t$r->[6], $pv\t$r->[5]\n"; 
      print "$s->[0]\t$s->[1]\t$r->[6], $pv\t$r->[5]\n";  
      
    } else {
      print OUT2 "$s->[0]\t$s->[1]\t\n";
      print "$s->[0]\t$s->[1]\t\n";
    }
    
    if (defined($outfile)) {
      foreach my $r (@$a_ref) {
	my $pv = sprintf("%1.2e", $r->[0]);
	print OUT "$s->[0]\t$s->[1]\t$pv\t$r->[1]\t$r->[2]\t$r->[3]\t$r->[4]\t$r->[6]\n";
      }
    }
    
  }
}

close OUT if (defined($outfile));
close OUT2 if (defined($outgo));




