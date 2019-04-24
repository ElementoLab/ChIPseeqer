use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


foreach my $r (@$a_ref) {

    print join("\t", @$r); print "\n";
    $r->[1] =~ s/\s//g;
    $r->[2] =~ s/\s//g;
    $r->[3] =~ s/\s//g;

    my $s1 = Sets::readSet("$ARGV[2]/CONSERVED_SETS/$r->[1].txt");
    my $s2 = Sets::readSet("$ARGV[2]/CONSERVED_SETS/$r->[2].txt");

    my $n1 = scalar(@$s1);
    my $n2 = scalar(@$s2);
    
    
    
    system("hypergeom -i $r->[3] -s1 $n1 -s2 $n2 -N $ARGV[1]");
    
} 
