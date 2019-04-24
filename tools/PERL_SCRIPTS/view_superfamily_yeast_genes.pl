#
#   display all the yeast genes 
#

use lib qw(/home/olly/PERL_MODULES);
use Superfamily;
use Database;


# get all yeast genes associated with the given SFAMID


my $db = Database->new;
$db->connect("YEASTHOUSE");

#$s = "select * from YEASTHOUSE where ORF='$ARGV[0]' or GENENAME='$ARGV[0]'";


#my $a = $db->queryOneRecord($s);
    
#$db->disconnect;

#print $a->{ORF} . "\n";

my $s = Superfamily->new;

my $a_ref = $s->getAllGenesBySFAMID($ARGV[0], 'sc');

foreach my $d (@$a_ref) {
    #print $d->{SEQUENCEID};

    #print "\t";
    my $s = "select * from YEASTHOUSE where ORF='$d->{SEQUENCEID}'";

    
    my $a = $db->queryOneRecord($s);
    
    #$db->disconnect;

    $d->{GENENAME} = $a->{GENENAME};

    #print "\n";
}

my @aa = sort { $a->{GENENAME} cmp $b->{GENENAME} } @$a_ref;


foreach my $d (@aa) {
    print $d->{SEQUENCEID};

    print "\t";
    
    print $d->{GENENAME};

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
