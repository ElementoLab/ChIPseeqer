BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $re = $ARGV[1];
my $rc = Sets::getREComplement($re);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    $s =~ s/$re/NNNNNNN/g;
    if (!defined($ARGV[2])) {
      $s =~ s/$rc/NNNNNNN/g;
    }
    print ">$n\n$s\n\n";
}

