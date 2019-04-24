use lib qw(/home/elemento/PERL_MODULES);

use  ClustalW;
use Sets;


my $cl = ClustalW->new;
$cl->setFile($ARGV[0]);

my $a_ref_aln = $cl->getSeqArray;

my @I = ();
my $l = undef;
my $n = scalar(@$a_ref_aln);
foreach my $seq (@$a_ref_aln) {
    my @a = split //, $seq;
    push @I, \@a;
    $l    = length($seq);
}


for (my $i=0; $i<$l; $i++) {
    my %H = ();
    for (my $j=0; $j<$n; $j++) {
	$H{ $I[$j][$i] } = 1;
    }
    
    my @a = keys(%H);
    if ((scalar(@a) == 1) &&
	(!defined($H{"-"}))) {
	print $a[0];
    } else {
	print "-";
    }
}
print "\n";
