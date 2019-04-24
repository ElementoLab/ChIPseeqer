#!/usr/bin/perl
#
#  usage : prg thr 7mers 8mers 9mers ..
#
use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;
use strict;

my $t = Table->new;

#   load the top 1000 7mers
#$t->setLimit(1000);

$t->loadFile($ARGV[0]);
my $a_refx = $t->getArray;


$t->setLimit(1000);

$t->loadFile($ARGV[1]);
my $a_refy = $t->getArray();


# get only kmers above last kmer's score
my $lowest = $a_refx->[ scalar(@$a_refx)-1 ]->[4];
print "lowest = $lowest\n";
my $a_refz = [];
foreach my $r (@$a_refy) {
    if ($r->[4] > $lowest) {
	push @$a_refz, $r;
    }
}



# for each of the new kmers, we want to make sure
# that there are no substrings in the old set
foreach my $r (@$a_refz) {
    
    # is there a substring in the refx list
    my $yes = undef;
    foreach my $s (@$a_refx) {
	
	#print "is $s->[0] a substring of $r->[0] ?\n";
	if (Sets::isSubstringOf_2S($s->[0], $r->[0]) || Sets::isSubstringOf_2S($r->[0], $s->[0])) {
	    $yes = $s->[0];
	}
    }

    

    if (!defined($yes)) {
	print "NO SUBS .." . join("\t", @$r) . "\n";
    } else {
	print "$r->[0] has a substring $yes ..\n";
    }
}

