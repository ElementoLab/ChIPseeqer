use lib qw(/home/olly/PERL_MODULES);
use Sets;


#Sets::printSet(Sets::allkmers($ARGV[0]));

#Sets::printSet(

my $a_ref = Sets::removeComplements(Sets::allkmers($ARGV[0]));

foreach my $r (@$a_ref) {
    my $H = Sets::countNucleotides($r);

    my @a = values(%$H);
    
    #print "a=" . join("\"", @a) . "\n";
    
    #print Sets::max(\@a);

    if (Sets::maxInArray(\@a) >= ($ARGV[0]-2)) {
	print "$r\n";
    }
}

#Sets::printSet(Sets::removeComplements(Sets::allkmers($ARGV[0])));


