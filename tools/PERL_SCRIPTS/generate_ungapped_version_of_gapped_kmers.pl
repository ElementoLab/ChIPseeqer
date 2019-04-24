use lib qw(/home/olly/PERL_MODULES);
use Sets;

$a_ref->[2] = Sets::allkmers(2);
$a_ref->[3] = Sets::allkmers(3);
$a_ref->[4] = Sets::allkmers(4);


my @BIN = ();

my %H = ();

open IN, $ARGV[0];
while (my $l = <IN>) {
    my @a = split /\t/, $l, -1;

    next if (length($a[0]) != $ARGV[1]);
    
    my ($g) = $a[0] =~ /(N+)/;
    my $k   = length($g);

    foreach my $r (@{ $a_ref->[$k] }) {
	my $n = $a[0];
	$n =~ s/N+/$r/;
	
	if (!defined($H{$n})) {
	    print "$n\n";
	    $H{$n} = 1;
	}
    }

	     

}


close IN;
