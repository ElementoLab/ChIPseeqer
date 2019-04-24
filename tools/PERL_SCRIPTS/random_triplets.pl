use lib qw(/home/olly/PERL_MODULES);
use Sets;
srand;
my $a_ref = Sets::readSet($ARGV[1]);

for (my $i=0; $i<$ARGV[0]; $i++) {
    
    my $r = Sets::sampleSetNoRep(3, $a_ref);

    print join("\t", @$r); print "\n";

    

}
