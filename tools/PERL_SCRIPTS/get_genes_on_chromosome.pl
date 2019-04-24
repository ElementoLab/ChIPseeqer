#
#  takes as input a set of genes, returns only those on chr2
#

use lib qw(/home/olly/PERL_MODULES);
use Database;
use Sets;

my $a_ref = Sets::readSet($ARGV[0]);


my $db = Database->new();
$db->connect("FLYDB");


#  get gene info
my $set         = Sets::SetToSQLSet($a_ref, "'");
my $sql         = "select * from FLY where ID in $set and CHROMOSOME like '$ARGV[1]'";
#print "$sql\n";
my $a_ref_genes = $db->queryAllRecordsRef($sql);

foreach my $r (@$a_ref_genes) {
    print "$r->{ID}\n";
}
