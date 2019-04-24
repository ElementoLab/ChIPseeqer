#
#  input: list of kmers, conserved sets, mapping
#

use lib qw(/home/olly/PERL_MODULES);
use Sets;


use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[3]);
my $a_ref_1_2 = $ta->getIndexKV(1, 0);
my $a_ref_2_1 = $ta->getIndexKV(0, 1);




# usage : script.pl kmers dir1 dir2 orthologs

my $a_ref_kmers = Sets::readSet($ARGV[0]);

foreach my $k (@$a_ref_kmers) {

    my $f1   = "$ARGV[1]/$k.txt";
    my $f2   = "$ARGV[2]/$k.txt";
    
    next if ((! -e $f1) || (! -e $f2));

    my @set1_o = ();
    my @set2_o = ();

    my $set1 = Sets::readSet($f1);
    foreach my $r (@$set1) {
	if (defined($a_ref_1_2->{$r})) {
	    push @set1_o, $r;
	}
    }


    my $set2 = Sets::readSet($f2);

    foreach my $r (@$set2) {
	if (defined($a_ref_2_1->{$r})) {
	    push @set2_o, $a_ref_2_1->{$r};
	}
    }

    my $a_ref_ov = Sets::getOverlapSet(\@set1_o, \@set2_o);

    my $s1 = scalar(@set1_o);
    my $s2 = scalar(@set2_o);

    if (scalar(@$a_ref_ov) > 0) {
	print join("\n", @$a_ref_ov);  print "\n";
    }
    print "$k\t$s1\t$s2\t"; print scalar(@$a_ref_ov); print "\n";
    
    
}

