BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Sets;
use strict;


#
# load up k-mers
#
my $a_ref_mo = Sets::readSet($ARGV[0]);

#
# get motif name
#
my $name = $ARGV[1];


#
# traverse sequences
#
my $fa = Fasta->new;
$fa->setFile($ARGV[2]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    #print ">$n\n$s\n\n";

    my @a_pos = ();

    foreach my $m (@$a_ref_mo) {
      my $a_ref = Sets::getREMotifPositionsOrientations($m, $s);
      push @a_pos, @$a_ref;
    }	

    
    foreach my $o (@a_pos) {
      print "$n\t$name\t$o->[0]\t$o->[1]\n";
    }
    
}
