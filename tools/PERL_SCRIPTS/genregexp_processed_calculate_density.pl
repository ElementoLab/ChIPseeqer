BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Fasta;
use Table;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my %L = (); 
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  $L{ $n } = length( $s );
}




my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

my %M = ();
foreach my $r (@$a_ref) {
  my $t = shift @$r;
  my $n = shift @$r;
  print "$t\t" . sprintf("%3.2f\n", 1000*$n/$L{$t});
}


