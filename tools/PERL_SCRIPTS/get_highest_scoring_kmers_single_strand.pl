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
$t->setLimit(1000);

my @aa = ();
for (my $i=1; $i<scalar(@ARGV); $i++) {
    
    # get the 1000 top kmers
    $t->loadFile($ARGV[$i]);
    my $a_refx = $t->getArray;
    
    
    # only add the kmers above threshold
    foreach my $r (@$a_refx) {
	push @aa, $r if (($r->[4] > $ARGV[0]) || $r->[4] =~ /inf/); 
    }
    
}



#   sort the kmers
my @bb = sort { $a->[4] <=> $b->[4] } @aa;

#   put the inf on top
while ($bb[0]->[4] =~ /inf/) {
    my $rr = shift @bb;
    push @bb, $rr;
}

# for each kmer, remove all substrings
my @bb = reverse(@bb);
my $n  = scalar(@bb);

for (my $i=0; $i<$n; $i++) {
    #print ">$bb[$i]->[0]\n";
    for (my $j=$i+1; $j<$n; $j++) {
	my $ssj = $bb[$j]->[0];


	if ((length($bb[$i]->[0]) == length($bb[$j]->[0])+1) &&
	   ($bb[$i]->[0] =~ /$ssj/)) {
	    $bb[$j]->[5] = "N";
	}

    }

}


foreach my $r (@bb) {
    if ($r->[5] ne "N") {
	print join("\t", @$r); print "\n";
    }
}


