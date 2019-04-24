BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Fasta;
use Sets;


my $h_ref = Sets::getIndex($ARGV[0]);

my $fa    = Fasta->new;
$fa->setFile($ARGV[1]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;

  if (defined($h_ref->{$n})) {
    print ">$n\n$s\n\n";
  }
  
}
