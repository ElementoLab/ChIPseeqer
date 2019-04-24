use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;


my $ta = Table->new;
$ta->setLimit($ARGV[1]);
$ta->loadFile($ARGV[0]);


my $a_ref = $ta->getArray();
foreach my $r (@$a_ref) {
    
    $H{ $r->[0] } ++;
    $H{ $r->[1] } ++;
    

}


Sets::printHash(\%H, "\t");



