# simple sql query 


#!/usr/bin/perl
# takes as input new col, + STDIN ORF tab feature
use lib qw(/home/olly/PERL_MODULES);
use Database;;
use Getopt::Long;
use Sets;
use GO_categories;
use DataFiles;
use Hypergeom;
use Table;

my $ta = Table->new;
#my $a_ref = Sets::readSet($ARGV[0]);
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getColumn(0);
# table !

my $df      = DataFiles->new;

my $db = Database->new();
$db->connect("FLYDB");

my $ntotal = 13522;

#  get gene info
my $set         = Sets::SetToSQLSet($a_ref, "'");
my $sql         = "select * from FLY where ID in $set;";
my $a_ref_genes = $db->queryAllRecordsRef($sql);
my $h_ref       = Sets::SQLRefToSimpleHash($a_ref_genes, "ID", "NAME");

foreach my $r (@$a_ref) {
    print "$r\t$h_ref->{$r}\n";
}


