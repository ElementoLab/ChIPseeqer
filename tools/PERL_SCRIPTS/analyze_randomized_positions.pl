use lib qw(/home/olly/PERL_MODULES);

use Sets;
use strict;

my $a_ref_files  = Sets::readSet($ARGV[1]);
my $a_ref_counts = Sets::getIndexKV($ARGV[0], 0, 1);

my %H = ();

foreach my $f (@$a_ref_files) {
    
    
    my $a_ref_counts_rand = Sets::getIndexKV($f, 0, 1);

    # add counts for each m
    foreach my $k (keys(%$a_ref_counts)) {
	push @{ $H{$k} }, $a_ref_counts_rand->{$k};	
    } 


}

my $nr = scalar(@$a_ref_files);

#  for each m, calculate the number of random counts above/equal actual value

foreach my $k (keys(%$a_ref_counts)) {
    
    #  get the actual count
    my $ac = $a_ref_counts->{$k};

    #  count 
    my $cnt = 0;
    
    foreach my $rc (@ { $H{$k} }) {
	if ($rc >= $ac) {
	    $cnt ++;
	}
    }

    print "$ac\t$cnt/$nr\n";
    
}
