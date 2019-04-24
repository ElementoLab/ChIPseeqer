use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;


my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();

my $n     = scalar(@$a_ref);

my @convergent_genes = ();
for (my $i=0; $i<$n; $i++) {

   
    my $r  = $a_ref->[$i  ];
    my $rp = $a_ref->[$i-1];

    if ($r->[0] eq $rp->[0]) {
	
	# add to the group
	push @group, $i;

    } else {

	if (scalar(@group) > 0) {
	    

	    # create a list of events (S, E)
	    my @list = ();
	    foreach my $m (@group) {
		my @a_tmp = ($m, "S", Sets::min($a_ref->[$m]->[1], $a_ref->[$m]->[2]));
		
		push @list, \@a_tmp;
		
		my @a_tmp = ($m, "E", Sets::max($a_ref->[$m]->[1], $a_ref->[$m]->[2]));

		push @list, \@a_tmp;
		

	    }

	    
	    my @list =  sort { $a->[2] <=> $b->[2] } @list;
	    
	    my $nl = scalar(@list);
	    for (my $j=0; $j<$nl-1; $j++) {
		
		
		
		#  ok ONLY if it an END
		next if ($list[$j]->[1] ne "E");

		# find the first START
		my $k  = $j + 1;
		my $S  = undef;
		while ($k < $nl) {
		    
		    # stop if found a further E
		    if ( ($list[$k]->[1] eq "E") && 
			 ($list[$k]->[2] > $list[$j]->[2]) ) {
			last;
		    } 

		    
		    # more complex if encounter a S
		    if ($list[$k]->[1] eq "S") {

			if (!defined($S) || ($list[$k]->[2] == $S)) {
			    # new pair if found 
			    
			    #if (($a_ref->[ $list[$j]->[0] ]->[3] == 1) && 
			    #	($a_ref->[ $list[$k]->[0] ]->[3] == -1)) {
				
			    if ($a_ref->[ $list[$j]->[0] ]->[3] == $a_ref->[ $list[$k]->[0] ]->[3]) {

				push @convergent_genes, $a_ref->[ $list[$j]->[0] ]->[4] if ($a_ref->[ $list[$j]->[0] ]->[3] ==  1);
				push @convergent_genes, $a_ref->[ $list[$k]->[0] ]->[4] if ($a_ref->[ $list[$k]->[0] ]->[3] == -1);

				$S = $list[$k]->[2];
			    }

			} else {
			    last;
			}
			
		    } else {
			last;
		    }


		    $k++;
		}

		# if the next token is a farther end, stop
		#next if (($list[$j+1]->[1] eq "E") && 
		#	 ($list[$j+1]->[2] > $list[$j  ]->[2]));

		# include all the next 
		
		#for (my $k=$j+1; $k<scalar(@list); $k++) {
		    
		    
		    
		#}
	    }

	}

	@group = ();

	push @group, $i;
    }

}

Sets::printSet(\@convergent_genes);
    
