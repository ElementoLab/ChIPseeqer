use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;

srand;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();


my $r = shift @$a_ref;
print join("\t", @$r); print "\n";

my $a_ref_shu = Sets::shuffle_array($a_ref);

foreach my $r (@$a_ref_shu) {
    
    

    print join("\t", @$r); print "\n";
    
    
    
}


