BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Sets;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
 
    
    print "\nAligning $n to genome ... \n";
    
    my $tmpfile = Sets::getTempFile("/tmp/titi");

    open OUT, ">$tmpfile";
    print OUT ">$n\n$s\n";
    close OUT;

    my $todo = "SIM4/sim4.2002-03-03/sim4.exe $tmpfile $ARGV[1]";
    
    system("$todo");

    
    
   
}
