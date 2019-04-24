BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use strict;

my $fa = Fasta->new;


$fa->setFile($ARGV[1]);
my %H1 = ();

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    $H1{ $n } = $s;
}
$fa->dispose();

open OUT1, ">f1.fa";
open OUT2, ">f2.fa";

my $fa = Fasta->new;

$fa->setFile($ARGV[0]);
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    #print "$n\n";

    if (defined( $H1{ $n } )) {
      print OUT1 ">$n\n$s\n\n";
      print OUT2 ">$n\n$H1{$n}\n\n";
    }
}

close OUT1;
close OUT2;
