BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use strict;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

open OUT, ">$ARGV[0].dict";

my %H = ();
my $cnt = 10000;
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  print ">$cnt\n$s\n\n";

  print OUT "$cnt\t$n\n";

  $cnt++;
}

close OUT;
