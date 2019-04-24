use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;



#
#  load transposons
#
my %TYPE = ();
open IN, $ARGV[1];
while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l;
    push @{ $TYPE{ $a[0] } }, \@a;
}
close IN;

my $a_ref_ty = Sets::readSet("transposon_types.txt");

print "\t" . join("\t", @$a_ref_ty); print "\n"; 

#
#  load genes
#
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
    print "$r->[0]";
    foreach my $type_t (@$a_ref_ty) {
	
	my $cnt = 0;
	foreach my $t (@{ $TYPE{ $type_t } }) {
	    
	    if (Sets::sequencesOverlap($r->[1]-100000, $r->[2]+100000, $t->[1], $t->[2])) {
		$cnt ++;
	    }
	}
	
	print "\t$cnt";

    }
    print "\n";

}
