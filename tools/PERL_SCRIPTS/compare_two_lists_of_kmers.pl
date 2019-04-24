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
my $a_ref = $t->getArray();

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


my $c2 = 0;
my $c3 = 0;
my $c4 = 0;



foreach my $r (@$a_ref) {

    print "1:" . join("\t", @$r);
    print "\t";
    if ($h_ref1->{ $r->[0] }) {
	$c2++;
	print "2:" . join("\t", @{ $h_ref1->{ $r->[0] } });
    } else {
	print "\t" x 1;
    }

    if ($h_ref2) {
	print "\t";
	if ($h_ref2->{ $r->[0] }) {
	    $c3++;
	    print "3:" . join("\t", @{ $h_ref2->{ $r->[0] } });
	} else {
	    print "\t" x 1;
	}
    }

    if ($h_ref3) {
	print "\t";
	if ($h_ref3->{ $r->[0] }) {
	    $c4++;
	    print "4:" . join("\t", @{ $h_ref3->{ $r->[0] } });
	} else {
	    print "\t" x 1;
	}
    }

    print "\n";
    
    
}

print "c2=$c2\nc2=$c3\nc4=$c4\n";
