use lib qw(/home/elemento/PERL_MODULES);
use Fasta;
use Sets;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my %H = ();
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
$s = uc($s);

    $s =~ s/U/T/g;

    $s = uc($s);

    my $k0 = substr($s, 0, 7);
    my $k1 = substr($s, 1, 7);
    #my $k2 = substr($s, 2, 7);

    push @{ $H{ Sets::getComplement(uc($k0)) } }, $n;
    push @{ $H{ Sets::getComplement(uc($k1)) } }, $n;
    #push @{ $H{ Sets::getComplement(uc($k2)) } }, $n;
    
    #print "$k0\t$k1\t$k2\n";
    

}


foreach my $k (keys(%H)) {
    print "$k\t";
    print join("\t", @{ $H{$k }}); 
    print "\n";
}
