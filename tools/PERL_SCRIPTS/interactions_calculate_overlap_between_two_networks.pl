use lib qw(/home/elemento/PERL_MODULES);

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref1 = $ta->getArray();

my %INDEX = ();
foreach my $r (@$a_ref1) {
  $INDEX{ $r->[0] }{ $r->[1] } = 1;
  $INDEX{ $r->[1] }{ $r->[0] } = 1;  
}

$ta->loadFile($ARGV[1]);
my $a_ref2 = $ta->getArray();

foreach my $r (@$a_ref2) {
  if ( defined($INDEX{ $r->[0] }{ $r->[1] }) ) {
    print "$r->[0]\t$r->[1]\n";
  }
}
