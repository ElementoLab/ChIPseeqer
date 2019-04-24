BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Sets;
use strict;


my $a_ref_files = Sets::getFiles($ARGV[0]);

foreach my $r (@$a_ref_files) {
  
  my @a = split /[\_\.]/, $r;
  
  

  my $fa = Fasta->new;
  $fa->setFile($r);
  
  while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    print "$n\t$a[2]\t$a[0]_$a[1]\n";
    
  }
}
