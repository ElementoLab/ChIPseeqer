BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Sets;
use strict;

my $fa = Fasta->new;
$fa->setFile($ARGV[1]);

my $a_ref_kmers = Sets::readSet($ARGV[0]);


while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  $s = lc($s);

  foreach my $r (@$a_ref_kmers) {
    my $se              =  $r;
    my $suc             =  uc($se);
    $s =~ s/$suc/$se/ig;
    $se              =  Sets::getComplement($se);
    $suc             =  uc($se);
    $s =~ s/$suc/$se/ig;
  }
  
  print ">$n\n$s\n\n";

}
