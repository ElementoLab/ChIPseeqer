use lib qw(/home/elemento/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %H = ();
foreach my $r (@$a_ref) {
  #$r->[3] = 1;
  if (!defined($H{ $r->[0] }{ $r->[1] }{ $r->[3] }) && !defined($H{ $r->[1] }{ $r->[0] }{ $r->[3] })) {
    print join("\t", @$r); print "\n";
    $H{ $r->[0] }{ $r->[1] }{ $r->[3] } = 1;
    $H{ $r->[1] }{ $r->[0] }{ $r->[3] } = 1;
  }


  
}
