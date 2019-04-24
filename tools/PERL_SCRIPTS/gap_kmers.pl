use lib qw(/home/olly/PERL_MODULES);
use Sets;

my $a_ref = Sets::readSet($ARGV[0]);

foreach my $r (@$a_ref) {
    
    my $ss1 = substr($r, 0, length($r)/2);
    my $ss2 = substr($r, length($r)/2, length($r)/2);

    print $ss1; print "." x $ARGV[1]; print "$ss2\n";
    
    
}
