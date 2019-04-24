use lib qw(/home/olly/PERL_MODULES);
use Sets;
use GroupEnrichment;
use strict;

my $go = GroupEnrichment->new;

$go->setGroups($ARGV[1]);

if ($ARGV[2] != -1) {
    $go->setGroupDesc($ARGV[2]);
}

my $a_ref_genes = Sets::readSet($ARGV[0]);
foreach my $r (@$a_ref_genes) {
    my $a = $go->getAnnotation($r);
    if ($a != -1) {
	print "$r\t"; print join("/", @$a); print "\n"; 
    }
}


