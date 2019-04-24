use lib "$ENV{FIREDIR}/SCRIPTS";

use Table;
use Sets;

srand;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();
my $r = shift @$a_ref;
print join("\t", @$r); print "\n";

my $a_ref_col = $ta->getColumn(0);
#shift @$a_ref_col;

my $a_ref_shu = Sets::shuffle_array($a_ref_col);

my $i = 0;
foreach my $r (@$a_ref) {
    
    $r->[0] = $a_ref_shu->[$i];

    print join("\t", @$r); print "\n";
    
    
    $i ++;
}


