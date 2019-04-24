use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;

#
#  input the top 500 gapped kmers
#
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


my $n = scalar(@$a_ref);

my @eliminate = ();

for (my $i=0; $i<$n-1; $i++) {
    
    for (my $j=$i+1; $j<$n; $j++) {
	
	next if (length($a_ref->[$i]->[0]) != length($a_ref->[$j]->[0]));

	my @a = split //, $a_ref->[$i]->[0];
	my @b = split //, $a_ref->[$j]->[0];
	my @c = ();
	my $l = length($a_ref->[$i]->[0]);
	my $out = 0;
	for (my $k=0; $k<$l; $k++) {
	    

	    if (($a[$k] eq 'N') || ($b[$k] eq 'N')) {
		$c[$k] = '.';
	    } elsif ($a[$k] eq $b[$k]) {
		$c[$k] = $a[$k];
	    } else {
		$out = 1; break;
	    }

	    
	}
	
	next if ($out == 1);
	
	my $pat = join("", @c);

	#print "$pat\n";
	
	if (($a_ref->[$i]->[0] =~ /$pat/) && ($a_ref->[$j]->[0] =~ /$pat/)) {
	    $eliminate[$j] = 1;
	}

    }

}

for (my $i=0; $i<$n; $i++) {
    #print "$eliminate[$i]\n";
    if (!defined($eliminate[$i])) {
	print join("\t", @{ $a_ref->[$i] }); print "\n";
    }
}
