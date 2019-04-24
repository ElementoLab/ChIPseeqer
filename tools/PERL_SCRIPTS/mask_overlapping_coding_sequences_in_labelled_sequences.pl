use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Fasta;
use strict;

my %CHR          = ();
my %BOUNDARIES_P = ();
my %STRAND       = ();
my %BOUNDARIES_T = ();

Sets::loadGeneTable($ARGV[1], \%CHR, \%BOUNDARIES_P, \%STRAND, \%BOUNDARIES_T);


my @keys = keys(%CHR);


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;

    my @a = split / /, $n, -1;

    my @a_exons = ();

    #if ($a[4] == -1) {
#	$s = Sets::getComplement($s);
#    }

    # get the overlapping transcript
    foreach my $k2 (@keys) {
	#next if ($a[0] eq $k2);

	#print "$a[0] ne $CHR{$k2}\n";
	
	next if (lc($a[0]) ne lc($CHR{$k2}));

	my @b = ($BOUNDARIES_T{ $k2 }->[0], $BOUNDARIES_T{ $k2 }->[1]);

	push @a_exons, \@b;
	
	#print "masking " . $BOUNDARIES_T{$k2}->[0] . ", " . $BOUNDARIES_T{ $k2 }->[1] . "\n";

	#next if (!Sets::sequencesOverlap($BOUNDARIES_T{ $a[0] }->[0], $BOUNDARIES_T{ $a[0] }->[1], $BOUNDARIES_T{ $k2 }->[0], $BOUNDARIES_T{ $k2 }->[1]) );

	# get the coordinate of the 
	#my $start_ovl = Sets::max(0, $BOUNDARIES_T{ $k2 }->[0] - $BOUNDARIES_T{ $a[0] }->[0]);
	#my $end_ovl   = Sets::min(length($s), $BOUNDARIES_T{ $k2 }->[0] -  $BOUNDARIES_T{ $a[0] }->[0]);

	
    } 

    my $s_masked = Sets::maskExons($s, \@a_exons, 'X');
#    if ($a[4] == -1) {
#	$s = Sets::getComplement($s);
#    }

    print ">$n\n$s_masked\n\n";
    
}
