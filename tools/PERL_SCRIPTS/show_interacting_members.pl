use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;

# input table of gene/miRNA interactions




my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndex(0);


# input table of interactions
open IN, $ARGV[1];

#$ta->loadFile($ARGV[1]);
#my $a_ref = $ta->getArray();

while (my $l = <IN>) {
    chomp $l; 
    my @c = split /\t/, $l;
    #foreach my $r (@$a_ref) {
    my $r = \@c;

    next if ($r->[0] eq $r->[1]);
    
    my @a = @{ $h_ref->{ $r->[0] } };
    my @b = @{ $h_ref->{ $r->[1] } };

    my $s = Sets::getOverlapSet(\@a, \@b); 
    
    my $n = scalar(@$s); 

    if ($n > 0) {
	print "$r->[0]\t$r->[1]\t$n\t"; 
	print join("\t", @$s); print "\n";
    }
}

close IN;




