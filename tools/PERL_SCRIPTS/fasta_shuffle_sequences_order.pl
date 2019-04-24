BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Sets;

srand(0);

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my %H = ();
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    $H{ $n } = $s;

}


my @kk = keys(%H);
my $kk_shu = Sets::shuffle_array(\@kk);

foreach my $k (@$kk_shu) {
  print ">$k\n$H{$k}\n\n";
}
