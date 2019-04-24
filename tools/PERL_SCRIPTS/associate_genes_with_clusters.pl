use lib qw(/home/olly/PERL_MODULES);
use Database;

my $db = Database->new();
$db->connect("FLYDB");


while (my $l = <STDIN>) {
    chomp $l;
    
    my @a = split /\t/, $l, -1;

    my $first = $a[1];
    my $last  = $a[$#a];
    
    my $sql = "select *  from FLY where CHROMOSOME = '$a[0]' AND ABS(ABS(START+END)/2 -ABS($first+$last)/2) < 10000";
	
    my $a_ref_genes = $db->queryAllRecordsRef($sql);

    foreach my $r (@$a_ref_genes) {
	print "$r->{ID}\t$r->{NAME}\n";
    }
}
