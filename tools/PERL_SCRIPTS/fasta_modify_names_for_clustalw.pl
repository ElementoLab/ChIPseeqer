BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);


my %H = ();
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  my $nn = substr($n, 0, 25);

  if (defined($H{ $nn })) {
    my $t = $H{ $nn };
    print ">$t" . "_$n\n$s\n\n";
    $H{ $nn } ++;
  } else {
    $H{ $nn } = 1;
    print ">$n\n$s\n\n";
  }
  
}
