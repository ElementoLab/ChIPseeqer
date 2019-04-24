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
  push @{ $M{ $r->[0] } }, $r->[2] if (!Sets::in_array($r->[2], @{ $M{ $r->[0] } }));
}

foreach my $k (keys(%M)) {
  my $n = scalar(@{ $M{ $k } });
  my $d = sprintf("%3.2f", 1000 * $n / $L{ $k }); 
  print "$k\t$d";  print "\n";
}

