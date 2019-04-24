use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;

my $f1    = shift @ARGV;
my $a_ref = Sets::readSet($f1);
my $ta    = Table->new;

my %H = ();
foreach my $f (@ARGV) {
    
    $ta->loadFile($f);
    my $a_ref1 = $ta->getArray();
    foreach my $r (@$a_ref1) {
	my @a = @$r; #print "$r->[0]\n";
	$H { $r->[0] } = \@a;
    }
    
    
}


foreach my $r (@$a_ref) {
    #print "kmer $r\n";
    print join("\t", @{ $H{ $r } });  print "\n";
    
}
