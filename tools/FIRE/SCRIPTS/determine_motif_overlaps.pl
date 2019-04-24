BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use Fire;

use strict;

# load summaryfile

# load list

# load profiles

my $tmpfile0 = "tmp-combined-motifs.txt";
Sets::writeSet(\@MOTIFS, $tmpfile0);

# if fastafile defined, create .profiles file ...
if (defined($fastafile)) {

  print "Generating profiles .. \n";
  my $todo = "perl $ENV{FIREDIR}/SCRIPTS/generate_motif_profiles.pl -fastafile $fastafile -summaryfile $tmpfile0 -rna $rna -noflank 1 -outfile $tmpfile0.profiles ";
  system("$todo");

  

}



my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %POS = ();
foreach my $r (@$a_ref) {
  push @{ $POS{ $r->[0] }->{ $r->[1] } }, $r->[2];
}


my @motifs = keys(%POS);

for (my $i=0; $i<@motifs-1; $i++) {

  my $tmp1       = Sets::get_array_from_re($motifs[$i]);
  my $k1         = @$tmp1;
  my $h_ref_gp1 = $POS{$motifs[$i]};
  my @g1        = keys( %$h_ref_gp1 );
  my $n1        = 0;

  foreach my $h (@g1) {
    $n1 += @{$h_ref_gp1->{$h}};      
  } 
  

  #print "num genes in g1 $h_ref_gp1 = " . scalar(@g1) . "\n";

  for (my $j=$i+1; $j<@motifs; $j++) {
    
    my $n2        = 0;
   
   # print "OVERLAP bet \n";

    my $tmp2       = Sets::get_array_from_re($motifs[$j]);
    my $k2         = @$tmp2;

    my $h_ref_gp2 = $POS{$motifs[$j]};
    my @g2        = keys( %$h_ref_gp2 );

    foreach my $h (@g2) {
      $n2 += @{$h_ref_gp2->{$h}};      
    } 
    
    my $gg        = Sets::getOverlapSet(\@g1, \@g2);

    #print "Overlap has " . scalar(@$gg) . "\n";

    my $ov        = 0;

    foreach my $g (@$gg) {
      
      foreach my $p1 (@{$h_ref_gp1->{$g}}) {
	
	foreach my $p2 (@{$h_ref_gp2->{$g}}) {
	  
	  if (Sets::sequencesOverlap($p1, $p1+$k1, $p2, $p2+$k2)) {
	    $ov ++;
	  }
	  
	}
	
      }
      
    }

    if (($ov/$n1>0.33) || ($ov/$n2>0.33)) {
      print "$motifs[$i] and $motifs[$j]: ov=$ov, n1=$n1, n2=$n2\n";
    }

  }
  
}



