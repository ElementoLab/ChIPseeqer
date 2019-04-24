use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;


my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref_kmers  = $ta->getArray();

my $n = scalar(@$a_ref_kmers);
for (my $i=0; $i<$n-1; $i++) { 

    my $k1 = $a_ref_kmers->[$i]->[0];
    
    print "-> $k1\n";

    for (my $j=$i+1; $j<$n; $j++) { 
	
	my $k2 = $a_ref_kmers->[$j]->[0];
    
	#print "$k1\t$k2\n";
	
	next if (length($k2)-length($k1) != 1); 

	if ($k2 =~ /^$k1/) {
	    #print "$k2\n"
	}

	if ($k2 =~ /$k1$/) {
	    print "$k2\n"
	}

    }

}

