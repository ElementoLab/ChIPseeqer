# suffle all sequences in a FASTA file
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;
use strict;

use Fasta;

my %H = ('A' => .24444, 'T' => .24444, 'C' => 0.25555, 'G' => 0.25555);




my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  my $news = "";
  for (my $i=0; $i<$ARGV[1]; $i++) {
    $news .= Sets::generateRandomSymbol(\%H); 
  }
  
  print ">$n\n$news\n\n";
  

}
