BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[1]);

my %H = ();
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    $H{$n} = $s;
}
$fa->dispose;

$fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {

  my ($n, $s) = @$a_ref;

  if (defined($H{$n})) {
    print ">$n\n$H{$n}\n\n";
  } else {
    print ">$n\n$s\n\n";
  }

}
