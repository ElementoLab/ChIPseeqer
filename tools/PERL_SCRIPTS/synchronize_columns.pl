use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;


my $ta = Table->new;
$ta->loadFile($ARGV[0]) or "Cannot load matrix1 file\n";
my $a_ref1 = $ta->getArray();

$ta->loadFile($ARGV[1]) or "Cannot load matrix2 file\n";
my $a_ref2 = $ta->getArray();

# get the overlap between the column IDs


my @a1 = @{ $a_ref1->[0] };
my @a2 = @{ $a_ref2->[0] };
shift @a1;
shift @a2;

my $a_ref_int = Sets::getOverlapSet(\@a1, \@a2);

#print "got " . scalar(@$a_ref_int) . " ov\n";

# matrix 1, what are the columsn to keep ?
my $i = 1;
my %H = ();
foreach my $c (@a1) {
    $H{ $c } = $i;
    $i ++;
}

my @tokeep = ("0");
foreach my $r (@$a_ref_int) {
    push @tokeep, $H{ $r };
}

#print "got " . scalar(@tokeep) . " cols to keeep\n";

open OUT, ">mm1.txt";
foreach my $r (@$a_ref1) {
    my @row = ();
    my $n = scalar(@tokeep);
    for (my $i=0; $i<$n; $i++) {
	push @row, $r->[ $tokeep[$i] ];
    }
    print OUT join("\t", @row); print OUT "\n";
}
close OUT;


# matrix 2, what are the columsn to keep ?
my $i = 1;
my %H = ();
foreach my $c (@a2) {
    $H{ $c } = $i;
    $i ++;
}

my @tokeep = ("0");
foreach my $r (@$a_ref_int) {
    push @tokeep, $H{ $r };
}


open OUT, ">mm2.txt";
foreach my $r (@$a_ref2) {
    my @row = ();
    my $n = scalar(@tokeep);
    for (my $i=0; $i<$n; $i++) {
	push @row, $r->[ $tokeep[$i] ];
    }
    print OUT join("\t", @row); print OUT "\n";
}
close OUT;


