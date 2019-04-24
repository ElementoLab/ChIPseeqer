#
#  align the second list on the first one
#

use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;

my $t = Table->new;

my $offset = 0;
if (int($ARGV[0]) != 0) {
    $offset = 1;
    $t->setLimit($ARGV[0]);
}


$t->loadFile($ARGV[0+$offset]);
my $a_ref = $t->getArray;

$t->loadFile($ARGV[1+$offset]);
my $h_ref1 = $t->getHash(0);

my $h_ref2 = undef;
if ($ARGV[2+$offset]) {
    $t->loadFile($ARGV[2+$offset]);
    $h_ref2 = $t->getHash(0);
}

my $h_ref3 = undef;
if ($ARGV[3+$offset]) {
    $t->loadFile($ARGV[3+$offset]);
    $h_ref3 = $t->getHash(0);
}



#
# go thru the 1rst reference
#
foreach my $r (@$a_ref) {
    
    my $out = 1;
    
    
    # is kmer in second file ?
    if (!defined($h_ref1->{ $r->[0]})) {
	next;
    }

    if ($h_ref2) {
	if (!defined($h_ref2->{ $r->[0] })) {
	    next;
	}
    }

    if ($h_ref3) {
	if (!defined($h_ref3->{ $r->[0] })) {
	    next;
	}
    }
    
    
    print $r->[0];
    print "\n";
    
    
}

