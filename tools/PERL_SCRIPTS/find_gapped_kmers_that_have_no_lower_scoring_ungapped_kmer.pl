use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;
use strict;

#
#  input the top 500 gapped kmers
#
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

#
#  input the gapped kmers here
#
my $a_ref_files = Sets::getFiles($ARGV[1]);

my $t           = $ARGV[2];

my $i = 1;
foreach my $r (@$a_ref) {

    #print "$r->[0] ?\n";

    my $fra_g = $r->[3] / $r->[1]; 

    my $len = length($r->[0]);
    
    my $pat = $r->[0];  $pat =~ s/N/\./g;

    #print "$r->[0]\t$r->[4]\n";

    # what is the best k-mer that matches this gapped one 
    my $bestscore = -10;
    my $bestkmer  = undef;
    my $bestfra  = undef;
    foreach my $f (@$a_ref_files) {
	next if (($f !~ /(\d+)mers_all\.txt/) && ($f !~ /(\d+)mers_ungapped_scored\.txt/));
	my ($k) = $1;

	next if ($k != $len);

	open IN, $f;
	#print "opening $f\n";
	while (my $l = <IN>) {
	    chomp $l;
	    my @a = split /\t/, $l, -1;

	    my $fra_u = $a[3] / ($a[1] + 0.0001); 

	    #  test if this kmer matched and scores higher than the gapped one
	    if (($a[0] =~ /$pat/) && ($a[4] > $bestscore)) {
		$bestscore = $a[4];
		$bestkmer = $a[0];
		$bestfra = $fra_u;
		#print "$i:" . $r->[0] . " -> " . $a[0] . "\t" . $r->[4] . " -> " . $a[4]; print "\n";
		
		#print "found $a[0], score = $bestscore\n";
	    } 

	    # if the current best score is higher than the score fror the gapped k-mer, we can skip
	    if ($bestscore > $r->[4]) {

		last;
	    }
	}

	close IN;

	if ($bestscore > $r->[4]) {
	    last;
	}
    
    }

    if (($bestscore < ($r->[4] - $t)) && ($fra_g > $bestfra)) {
	print join("\t", @$r); print "\t$bestkmer\t$bestscore\t(" . sprintf("%3.2f, %3,2f", $fra_g, $bestfra) . ")\n";
    }
    
    $i++;
}
