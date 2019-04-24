use lib qw(/home/olly/PERL_MODULES);

require "Mips.pm";

$m1 = Mips->new;

$m1->readORFset($ARGV[0]);


$m1->setPvalueThreshold(0.01);


$m2 = Mips->new;

$m2->readORFset($ARGV[1]);


foreach my $r (@{$m1->getFunctionalContent}) {
    print "$r->{NUM}\t$r->{PVALUE}\t$r->{OVERLAP} / $r->{TOTAL}\t$r->{TEXT}";
   
    $r2 = $m2->getSpecificEnrichment($r->{NUM});

    
    print "\t$r2->{PVALUE}\t$r2->{OVERLAP} / $r2->{TOTAL}";
    
    

    print "\n";
}


