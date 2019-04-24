use lib qw(/home/olly/PERL_MODULES);

use Sets;


my $a_ref = Sets::getArray($ARGV[0]);

my @a = ();
foreach my $r (@$a_ref) {
    push @a, $r->[0];
    push @a, $r->[1];
}


my $n = scalar(@$a_ref);
srand;
foreach my $r (@$a_ref) {
    my $idx1 = rand($n);
    my $idx2 = rand($n);
    
    print "$a[$idx1]\t$a[$idx2]\n";
}

