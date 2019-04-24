use lib qw(/home/olly/PERL_MODULES);

use Sets;
use Table;
use strict;
my $ta = Table->new;

my $a_ref_files = Sets::getFiles($ARGV[0]);

my @A = ();
foreach my $f (@$a_ref_files) {
    
    #print "$f\n";
    $ta->loadFile($f);
    
    my $a_ref = $ta->getColumn(4);

    my $per   = Sets::percentile($a_ref, 0.999);
    
    #my $avg = Sets::average($a_ref);
    
    #print "$per\n";

    push @A, $per;

}


print Sets::average(\@A); print "\n";
