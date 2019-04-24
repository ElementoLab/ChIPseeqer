use lib qw(/home/olly/PERL_MODULES /home/olly/PERL_MODULES/Hypergeom /home/olly/PERL_MODULES/Hypergeom/blib/lib /home/olly/PERL_MODULES/Hypergeom/blib/arch/auto/Hypergeom);
use GO_func;
use Sets;
use Table;

my %H = ();
my $ta = Table->new;
$ta->loadFile("/home/olly/DATA/YEASTS/SGD_features.tab");
my $a_ref = $ta->getArray();
foreach my $r (@$a_ref) {
    #print "$r->[1]\n";
    if ($r->[1] eq "CDS") {
	$H{ $r->[6] } = 1;
	#print "$r->[6]\n";
    }
}


open IN, $ARGV[0];
my $i = 0;
while (my $l = <IN>) {
    chomp $l;
    
    my @a = split /\t/, $l;

    if (defined($H{$a[0]}) && defined($H{$a[1]}) && defined($H{$a[2]})) {
	print "$l\n";
    }
}


   
