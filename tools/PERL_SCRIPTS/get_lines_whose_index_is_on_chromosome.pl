#
#  takes as input a set of genes, returns only those on chr2
#

use lib qw(/home/olly/PERL_MODULES);
use Database;
use Sets;
use Table;


my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


my $db = Database->new();
$db->connect("FLYDB");


#  get gene info
my $set         = Sets::SetToSQLSet($ta->getColumn(0), "'");
my $sql         = "select ID from FLY where ID in $set and CHROMOSOME like '$ARGV[1]'";

my $a_ref_genes = $db->queryAllRecordsRef($sql);
my $h_ref_chr   = Sets::SQLRefToIndex($a_ref_genes, "ID");

foreach my $r (@$a_ref) {
    
    if (defined($h_ref_chr->{$r->[0]})) {
	print join("\t", @$r); print "\n";
    }
}
