use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;
use strict;

my $a_ref_files = Sets::getFiles($ARGV[0]);

my @COUNTS = ();
foreach my $f (@$a_ref_files) {
    
    my $a_ref = Sets::readSet($f);

    foreach my $r (@$a_ref) {
	$COUNTS[ $r ] ++;
    }
}


my $n  = scalar(@COUNTS);
my $nf = scalar(@$a_ref_files);

for (my $i=0; $i<$n; $i++) {
    #my $t = int(0.5 + $COUNTS[$i] / $nf );

    my $t = $COUNTS[$i] / $nf;
    print "$i\t$t\n";
}



