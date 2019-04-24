BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;



my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndex(0);

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[1]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  #$n = "L$n";
  if ( defined($h_ref->{$n}) ) {
    print ">$n\n$s\n\n";
  }
}
