#
#  finds cluster from v2 files, then tries to link adjacent clusters .. (separated by X genes)
#


#
#  reads in the lists of duplicates (unique)
#

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;
use POSIX qw(ceil floor);


#
#  load the gene information (genes only, no proteins)
#
my $ta = Table->new;

#
# create the adjacency list
#
my %ADJ = ();

my %EDGE = ();

my $double = $ARGV[1];

#
#  then traverse the duplicates in search for those < d
#
$ta->setFile($ARGV[0]);
while (my $r = $ta->nextRow()) {

    if (defined($double)) {
	$EDGE{ $r->[0] } { $r->[1] } = 1;
	if (defined($EDGE{ $r->[1] } { $r->[0] })) {
	    push @{ $ADJ{ $r->[0] } }, $r->[1];
	    push @{ $ADJ{ $r->[1] } }, $r->[0];
	}
    } else {
	push @{ $ADJ{ $r->[0] } }, $r->[1];
	push @{ $ADJ{ $r->[1] } }, $r->[0];
    }
} 
$ta->dispose;


# init colors
my %COL = ();
foreach my $n (keys(%ADJ)) {
    $COL{$n} = 0;
}

my @size = ();
foreach my $n (keys(%ADJ)) {

    
    
    if ($COL{$n} == 0) {
	
	#print "Start at $n ?\n";
	
	my @a_nodes = ();
	DFS(\%ADJ, $n, \%COL, \@a_nodes);

	print join("\t", @a_nodes); print "\n";

	
	$size [ scalar(@a_nodes) ] ++;
	
	
	
    }
    
} 
	





	
	


my $nsize = scalar(@size);
for (my $i=0; $i<$nsize; $i++) {
    print "$i\t$size[$i]\n" if (defined($size[$i]));
    
}



sub DFS {
    my ($h_ref_adj, $node, $h_ref_col, $a_ref_nodes) = @_;

    $h_ref_col->{ $node } = 1;
    
    #print "Visiting $node ..\n";

    foreach my $next (@{ $h_ref_adj->{ $node } }) {

	next if (!defined($h_ref_adj->{$next}));


	if ($h_ref_col->{ $next } == 0) {
            DFS($h_ref_adj, $next, $h_ref_col, $a_ref_nodes);
        }
    }
    
    $h_ref_col->{ $node } = 2;

    push @$a_ref_nodes, $node;

}

