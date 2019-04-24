use lib qw(/home/elemento/PERL_MODULES);


use Fasta;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while ( my $a_ref = $fa->nextSeq ) {
    my ($name, $seq) = @{$a_ref};
    my @a = split /\ /, $name;
    if ($name =~ /$ARGV[1]/) {
	print ">$a[0]\n$seq\n";
    }
}
