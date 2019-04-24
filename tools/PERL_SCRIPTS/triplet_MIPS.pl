use lib qw(/home/olly/PERL_MODULES /home/olly/PERL_MODULES/Hypergeom /home/olly/PERL_MODULES/Hypergeom/blib/lib /home/olly/PERL_MODULES/Hypergeom/blib/arch/auto/Hypergeom);
use GO_func;
use Sets;

$m = GO_func->new;
$m->setSource("MIPS", "YEAST");
$m->setTotalNbORFS(5500);
#$m->setPvalueThreshold(0.00001);


$m->setORFset($a_ref);
$m->setVerbose(0);


open IN, $ARGV[0];
my $i = 0;
while (my $l = <IN>) {
    chomp $l;
    
    my @a = split /\t/, $l;

    my $r1 = $m->getMIPScategory($a[0]);
    my $r2 = $m->getMIPScategory($a[1]);
    my $r3 = $m->getMIPScategory($a[2]);

    my $i1 = Sets::getOverlapSet($r1, $r2);
    my $i2 = Sets::getOverlapSet($i1, $r3);
    
    if (scalar(@$i2) > 0) { 
	
	if ($i2->[0] !~ /UNCLASSIFIED PROTEINS/) {
	    print "$l\t"; print join("\t", @$i2); print "\n";
	}
    }

    $i ++;


    last if ($i == 10000);

}
close IN;
