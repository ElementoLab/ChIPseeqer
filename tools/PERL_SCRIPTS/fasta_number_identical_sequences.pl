BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use strict;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);


my %COUNTS = ();
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    my $nn = undef;
    if (defined($COUNTS{$n})) {
      $nn = "$n-$COUNTS{$n}"; 
    } else {
      $nn = $n;
    }

    $COUNTS{$n} ++;
    
    print ">$nn\n$s\n\n";
    
}
