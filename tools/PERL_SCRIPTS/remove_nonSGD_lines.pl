use lib qw(/home/olly/PERL_MODULES /home/olly/PERL_MODULES/Hypergeom /home/olly/PERL_MODULES/Hypergeom/blib/lib /home/olly/PERL_MODULES/Hypergeom/blib/arch/auto/Hypergeom);
use GO_func;
use Sets;
use Table;

my %H = ();
my $ta = Table->new;
$ta->loadFile("/home/olly/DATA/YEASTS/SGD_features.tab");
my $a_ref = $ta->getArray();
foreach my $r (@$a_ref) {
    if ($r->[1] eq "CDS") {
	$H{ $r->[6] } = 1;
    }
}


open IN, $ARGV[0];
my $i = 0;
my $l = <IN>;
print $l;
while (my $l = <IN>) {
    chomp $l;
    
    my @a = split /\t/, $l;

    if (defined($H{$a[0]})) {
	print "$l\n";
    }
}
close IN;

   
