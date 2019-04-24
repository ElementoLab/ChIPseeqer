BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Alignment;


my $al = Alignment->new;
$al->setComp('RNA');
$al->setVerbose(1);
my $a_ref = $al->nw("GCCTACCTCT", "TTGGTGTGTTGGATGATGGAGT");

#my $a_ref = $al->nw("TACGGT", "ATCAGCGTCA");
#my $a_ref = $al->nw("ATAAAAT", "ATAT");
#my $a_ref = $al->nw("MPRCLCQRJNCBA", "PBRCKCRNJCJA");

print "$a_ref->[0]\n";
print "$a_ref->[1]\n";

#$al->nw("ACGTTTACGGT", "AAGGTACGC");





#$al->nw("ATGATA", "ATCAGGGTCA");
#$al->nw("ATATATATA", "ATCAGTCA");



# -PWDNAMATRIX CLUSTALW -PWGAPOPEN=f -PWGAPEXT=f
