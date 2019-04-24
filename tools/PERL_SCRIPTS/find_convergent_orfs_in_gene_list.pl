use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;
use strict;


my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();

my $n     = scalar(@$a_ref);



my @convergent_genes = ();
my %CHRS = ();

my $chr_idx  = 1;
my $orf_idx  = 0;
my $sta_idx  = 2;
my $end_idx  = 3;
my $str_idx  = 4;

my $idx = 0;
foreach my $r (@$a_ref) {
    push @$r, $idx;
    push @{ $CHRS{ $r->[$chr_idx] } }, $r;
    $idx ++;
} 

foreach my $c (keys(%CHRS)) {
 
    
    #
    #  NEW CHROMOSOME
    #
    my @group = @{ $CHRS{ $c } };
    
    # create a list of events (S, E)
    my @list = ();
    foreach my $r (@group) {
	my @a_tmp = ($r->[$orf_idx], "S", Sets::min($r->[$sta_idx], $r->[$end_idx]), $r->[$str_idx]);
	push @list, \@a_tmp;
	my @a_tmp = ($r->[$orf_idx], "E", Sets::max($r->[$sta_idx], $r->[$end_idx]), $r->[$str_idx]);
	push @list, \@a_tmp;
    }

    # sort list of events
    my @list =  sort { $a->[2] <=> $b->[2] } @list;
    
    # go thru the sorted list, searching for convergent genes
    my $nl = scalar(@list);
    for (my $j=0; $j<$nl-1; $j++) {
		
	#  continues ONLY if it an END
	if (($list[$j  ]->[1] eq "E") && ($list[$j  ]->[3] ==  1) &&
	    ($list[$j+1]->[1] eq "S") && ($list[$j+1]->[3] == -1)) {


	    #next if (abs($list[$j  ]->[2] -  $list[$j+1]->[2]) < 1000);
	    
	    push @convergent_genes, $list[$j  ]->[0];
	    push @convergent_genes, $list[$j+1]->[0];

	    #print $list[$j  ]->[0] . " and " . $list[$j+1]->[0] . " are convergent ..\n";

	}
    }
}


Sets::printSet(\@convergent_genes);
    
