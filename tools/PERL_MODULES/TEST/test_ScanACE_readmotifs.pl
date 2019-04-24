BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use ScanACE;


my $sc = ScanACE->new;

#$sc->readAlignACEMotifs($ARGV[0]);
#my $a_ref = $sc->getMotifsTXT();
#foreach my $r (@$a_ref) {
#open OUT, ">m1.txt";
#close OUT;
#}



$sc->setMotif($ARGV[0]);
$sc->setFasta($ARGV[1]);
$sc->setStdDev(2.0);
$sc->setGC(0.44);
$sc->run;

my $cnt = $sc->getNbMotifsInFile($ARGV[0]);
for (my $i=0; $i<$cnt; $i++) {
  my $a_ref = $sc->getSites($i);
  foreach my $a (@$a_ref) {
    print "$ARGV[0]\t$a->[0]\t$a->[3]\t$a->[4]\t$a->[1]\t$a->[2]\n";
  }
}
