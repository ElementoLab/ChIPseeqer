use lib qw(/home/olly/PERL_MODULES);

use Sets;
use Table;


Sets::cmdLine(2, "ensembl2swissprot.pl GENES TRADUC");
#
#  read set of ensembl genes
#


#
#  read traduction ensembl -> swissprot
#


my $a_ref = Sets::readSet($ARGV[0]);


my $ta = Table->new;
$ta->loadFile($ARGV[1]);


my $h_ref_annot = $ta->getIndex(0);

my @a_set = ();
foreach my $r (@$a_ref) {
    
    if ($h_ref_annot->{$r}->[1]) {
	my @a = split /,/, $h_ref_annot->{$r}->[1];
	push @a_set, @a;
    }
    
}

my $a_ref_set_uniq = Sets::removeDuplicates(\@a_set);

Sets::printSet($a_ref_set_uniq);
