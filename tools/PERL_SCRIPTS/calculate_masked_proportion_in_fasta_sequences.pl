use lib qw(/home/olly/PERL_MODULES);
use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $cnt = 0;
my $HN  = 0;
my $HX  = 0;
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    my @a = split //, $s;
    
    foreach my $r (@a) {
        next if ($r eq 'N');
	$HN  ++ if ($r ne 'X');
	$HX  ++ if ($r eq 'X');
	$cnt ++;
    }

}

print "X : $HX, others : $HN, total : $cnt\n";
