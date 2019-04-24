BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
 
    my @a = split /\|/, $n;
    
    $a[4] =~ s/^\ +//g;
    $a[4] =~ s/\ \[Caulobacter\ crescentus\ CB15\]//g;
    print "$a[3]\t$a[4]\n";
   
}
