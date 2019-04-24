BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";




use  ClustalW;
use Sets;
use GenRegexp;
use DataFiles;



my $fasta = $ARGV[0];
my $a_ref_block = Sets::readSet($ARGV[1]);


foreach my $k (@$a_ref_block) {
    
    my $ge = GenRegexp->new;
    $ge->calc($k, $fasta);
    
    my $c = $ge->getNbMatches(); 
    
    print "$k\t$c\n";

}


