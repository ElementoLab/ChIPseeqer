use lib qw(/home/elemento/PERL_MODULES);

use Fasta;
use Sets;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $a_ref_aln = [];
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    push @{ $a_ref_aln }, uc($s);
}



my @I = ();
my $l = undef;
my $n = scalar(@$a_ref_aln);
foreach my $seq (@$a_ref_aln) {
    my @a = split //, $seq;
    push @I, \@a;
    $l    = length($seq);
}


if ($ARGV[1]) {
  print ">$ARGV[1]\n";
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
