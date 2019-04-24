#!/usr/bin/perl

use lib qw(/home/olly/PERL_MODULES);

use Sets;
use Table;

# load conserved set
my $a_ref_set = Sets::readSet($ARGV[0]);

# load llid / ensid table
my $t = Table->new;
$t->loadFile("/home/olly/DATA/HUMAN/ENSEMBL_GENES/ENSEMBLID_TO_LOCUSLINK/ensllidu.txt");

my $h_ref_idx = $t->getIndexColumnsKV(1, 0);

my @O = ();
foreach my $r (@$a_ref_set) {
        
    if (defined($h_ref_idx->{$r})) {
        push @O, $h_ref_idx->{$r};
    } 

    
    
}

my $a_ref_out = Sets::removeDuplicates(\@O);

Sets::printSet($a_ref_out);

