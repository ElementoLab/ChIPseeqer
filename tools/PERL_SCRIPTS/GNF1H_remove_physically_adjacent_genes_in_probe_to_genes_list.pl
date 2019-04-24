#
#  if a probe matches X genes that are in less than 10kb of each other, remove random one
#   how about > 2 genes ?
use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;
use GenomeAnnotation;
use strict;


my $ga = GenomeAnnotation->new;

$ga->read($ARGV[1]);


#print $ga->genomicDistanceBetweenGenes( "CG33105", "CG33104" );


#print "\n";



Sets::cmdLine(1, "probe2genes genes_infos");




my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %INDEX_REM = ();
my %H         = ();
foreach my $r (@$a_ref) {
    
    my @b = split /\//, $r->[1], -1;

    my $n = scalar(@b);

    if ($n > 1) {
	
	my $h_ref_ordered = $ga->orderGenesOnChromosomes(\@b);
	
	my @chrs          = keys(%$h_ref_ordered);

	#print "genes are on " . join(",", @chrs) . "\n";
	
	my @removed = ();
	foreach my $c (@chrs) {
	    
	    my $nc = scalar(@{ $h_ref_ordered->{$c} });

	    next if ($nc == 1);

	    #
	    # reverse the order if the first member if in direction -1
	    #
	    if ($ga->getStrand($h_ref_ordered->{$c}->[0]) == 1) {
		@{ $h_ref_ordered->{$c} } = reverse @{ $h_ref_ordered->{$c} };
	    }

	    for (my $i=0; $i<$nc-1; $i++) {
		
		# if the next gene is too 
		my $d = $ga->genomicDistanceBetweenGenes($h_ref_ordered->{$c}->[$i], $h_ref_ordered->{$c}->[$i+1]);
		
		my $t = $ga->getBoundariesCDS($h_ref_ordered->{$c}->[$i]); #print "s=$t->[0]\te=$t->[1]\n";

		#print "distance between " . $h_ref_ordered->{$c}->[$i] . " and " . $h_ref_ordered->{$c}->[$i+1] . " is $d\n";
		
		if (!defined($d)) {
		    die "problem\n";
		} elsif ($d < 10000) {
		    push @removed, $h_ref_ordered->{$c}->[$i];
		}
		
	    }
	}
	
	my @left = ();
	foreach my $g (@b) {
	    if (!Sets::in_array($g, @removed)) {
		push @left, $g;
	    } else {
		$INDEX_REM{ $g } = 1;
	    }
	}
	
	print "$r->[0]\t"; print join("\/", @left); print "\n";

	#<STDIN>;

    } else {
	print "$r->[0]\t"; print join("\/", @b); print "\n";
    }

}


open OUT, ">genes_removed.txt";
print OUT join("\n", keys(%INDEX_REM));
print OUT "\n";
close OUT;
