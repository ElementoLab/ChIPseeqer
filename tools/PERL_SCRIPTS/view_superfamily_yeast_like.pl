#
#   display all the yeast genes 
#

use lib qw(/home/olly/PERL_MODULES);
use Superfamily;
use Database;

my $db = Database->new;
$db->connect("YEASTHOUSE");

$s = "select * from YEASTHOUSE where ORF='$ARGV[0]' or GENENAME='$ARGV[0]'";


my $a = $db->queryOneRecord($s);
    
$db->disconnect;

print $a->{ORF} . "\n";

my $s = Superfamily->new;

my $a_ref = $s->getDomains($a->{ORF});

my @aa = sort  sortByRegion @$a_ref;

foreach my $d (@aa) {
    print ">>" . $a->{ORF}; print "\t";
    print $d->{REGION};
    print "\t";
    print $d->{SFAMID};
    print "\t";
    print $d->{TEXT};
    print "\n";
}


sub  sortByRegion {

    my $ar  = $a->{REGION};
    my $br  = $b->{REGION};

    my ($ars) = $ar =~ /^(\d+)\-/;
    my ($brs) = $br =~ /^(\d+)\-/;

    #print "$ars\t$brs\n";

    return $ars <=> $brs; 
}
