# 1. load lengths
use lib qw(/home/elemento/PERL_MODULES);
use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndexKV(0, 1);

# 2. go thru the list
my $a_ref = Sets::readSet($ARGV[0]);
my @a = ();
foreach my $r (@$a_ref) {
  push @a, $h_ref->{ $r } if (defined($h_ref->{ $r }));
}

print Sets::median(\@a); print "\n";

#Sets::printSet(\@a);
