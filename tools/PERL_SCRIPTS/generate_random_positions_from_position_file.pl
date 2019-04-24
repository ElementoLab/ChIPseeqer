use lib qw(/home/olly/PERL_MODULES);
use Sets;

srand;
my $a_ref = Sets::readSet($ARGV[0]);

my $max   = Sets::max($a_ref);
my $n     = scalar(@$a_ref);

my @a = ();
for (my $i=0; $i<$n; $i++) {
    
    push @a,  int(0.5+rand($max));
    #print "\n";
}



my @b = sort { $a <=> $b } @a;

Sets::printSet(\@b);
