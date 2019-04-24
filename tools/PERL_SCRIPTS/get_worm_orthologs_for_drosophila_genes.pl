#!/usr/bin/perl
#
# takes as input set of droso genes, mapping 
#
use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;

my $ta = Table->new;
$ta->loadFile("/home/olly/PROGRAMS/KMERS_ORTHOLOGS/MAPPINGS/ELE_MEL-final-drosofirst-fbgns.txt");


my $h_ref_map = $ta->getIndexKV(0, 1);

my $a_ref_dro = Sets::readSet($ARGV[0]);

foreach my $r (@$a_ref_dro) {
    if (defined($h_ref_map->{$r})) {
	print $h_ref_map->{$r} . "\n";
    }
}
