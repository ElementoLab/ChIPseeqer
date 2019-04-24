use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;
use Database;
use strict;
use Getopt::Long;

my $ids      = undef;
my $suffix   = undef;
my $database = undef;
my $table    = undef;
my $index    = undef;

GetOptions ('ids=s'      => \$ids,
	    'suffix=s'   => \$suffix,
	    'database=s' => \$database,
	    'table=s'    => \$table,
	    'index=s'    => \$index);

# get the set of genes
my $a_ref = Sets::readSet($ids);

my $db = Database->new();
$db->connect($database);


#
#  create a new table
#
my $sql = `mysqldump -d $database $table`;
my $bef = "CREATE TABLE $table";
my $aft = "CREATE TABLE $table" . "_$suffix";
$sql =~ s/$bef/$aft/;

my $file = Sets::saveToTempFile($sql);

print "mysql $database < $file\n";
system("mysql $database < $file");

unlink $file;

#
#  recreate data
#
my $tmpfile = Sets::getTempFile("/tmp/tci");
open OUT, ">$tmpfile";
foreach my $r (@$a_ref) {
    
    my $sql = "select * from $table where $index = '$r';";

    my $a_ref_rows = $db->queryAllRecordsArrayRef($sql);
    
    my $n = scalar(@$a_ref_rows);
    
    if ($n == 0) {
	print "No DB info for $r\n";
    }

    foreach my $s (@$a_ref_rows) {
	print OUT join("\t", @$s); print OUT "\n";
	#print join("\t", @$ss); print "\n";

	
    }
}
close OUT;
print "mysql $database -e \"load data local infile '$tmpfile' into table $table" . "_$suffix\"\n";
system("mysql $database -e \"load data local infile '$tmpfile' into table $table" . "_$suffix\"");

unlink $tmpfile;


