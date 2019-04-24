use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;


my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


my $a_ref_files = Sets::getFiles($ARGV[1]);

my $i = 1;
foreach my $r (@$a_ref) {

    my $len = length($r->[0]);

    #print "$r->[0]\t$r->[4]\n";

    foreach my $f (@$a_ref_files) {

	my ($k,$k1,$g) = $f =~ /^(\d)\.(\d)\.(\d)/;

	next if (($k+$g) != $len);

	open IN, $f;
	#print "opening $f\n";
	while (my $l = <IN>) {
	    chomp $l;
	    my @a = split /\t/, $l;

	    my $pat = $a[0];  $pat =~ s/N/\./g;

	    #print "-> $pat $a[4]\n"; 

	    

	    if (($r->[0] =~ /$pat/) && ($a[4] > $r->[4])) {
		print "$i:" . $r->[0] . " -> " . $a[0] . "\t" . $r->[4] . " -> " . $a[4]; print "\n";
	    } else {
		last;
	    }
	}

	close IN;
    }
    
    $i++;
}
