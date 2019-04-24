BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[1]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  my $p = length($s) - $ARGV[0];
  my $ss = substr($s, Sets::max(0, $p), $ARGV[0]);  
  
  print ">$n\n$ss\n\n";
}
