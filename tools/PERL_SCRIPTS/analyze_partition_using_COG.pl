use lib qw(/home/olly/PERL_MODULES);
use Sets;
use COG;
use Table;
use strict;

my $co = COG->new;

$co->setSpecies($ARGV[1]);
$co->setTotalNbORFS($ARGV[2]);
#$co->setBonferroni(1);

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %PARTITION = ();

foreach my $r (@$a_ref) {
    push @{ $PARTITION{  $r->[2] } },  $r->[1];

    #print "add $r->[0] to partition $r->[1]\n";
}

foreach my $p (keys(%PARTITION)) {
    
    print "COG enrichments:\n";
    #print join("\t", @{ $PARTITION{$p} });
    $co->setORFset($PARTITION{$p});
    my $a_ref = $co->getFunctionalContent();
    foreach my $c (@$a_ref) {
	print "\t$c->{TEXT}\t" . sprintf("%3.2e\n", $c->{PVALUE});
    }

}

