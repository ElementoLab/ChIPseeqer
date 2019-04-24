BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use ScanACE;
use Hypergeom;

use Fasta;
use Sets;

#
# 0. read list of interest
#
my $h_ref_e = Sets::getIndex($ARGV[1]);
my $s1      = scalar(keys(%$h_ref_e));

#
# 1. read enhancers
#
my $fa = Fasta->new;
$fa->setFile($ARGV[0]);
my $N = 0; 
my %L = ();

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  $L{ $n } = length( $s );
  $N ++;
}

#
# 2. scanace
#
my $sc = ScanACE->new;
$sc->setMotif($ARGV[2]);
$sc->setFasta($ARGV[0]);
$sc->setStdDev(3.0);
$sc->setGC(0.44);
$sc->run;

my $cnt = $sc->getNbMotifsInFile($ARGV[2]);
for (my $i=1; $i<=$cnt; $i++) {
  my $a_ref = $sc->getSites($i);
  $s2 = 0;
  $ov = 0;
  my %H = ();
  my $cnt = 0;

  #next if (length($a_ref->[0]->[3]) > 11);
  
  foreach my $a (@$a_ref) {
    print "$a->[3]\t" if ($cnt == 0);
    $cnt ++;
    $H{ $a->[0] }++;
  } 
  
  foreach my $k (keys(%H)) {
    if (($H{ $k } > 1) && (1000*$H{ $k }/$L{ $k } > 2.0)) {

      $s2 ++;
      if (defined($h_ref_e->{ $k })) {
	print "$k\n";
	$ov ++;
      }
    }
  }
  
    
    #if (!defined($H{ $a->[0] })) {
    #  $s2 ++;
    #  if (defined($h_ref_e->{ $a->[0] })) {
#	$ov ++;
#      }
#    }
#    $H{ $a->[0] } = 1;
#  }
  #print "Hypergeom::lcumhyper($ov, $s1, $s2, $N);\n";
  my $p = Hypergeom::cumhyper($ov, $s1, $s2, $N);    
  $p = sprintf("%4.3e", $p);
  print "$i\t$ov, $s1, $s2, $N\t$p\n";
}
