BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;

# read in operon table
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %H = ();
foreach my $r (@$a_ref) {
  my @a = split /\-/, $r->[0], -1;
  
  foreach my $s (@a) {
    $H{ $s } = $r->[0]; #print "$s => $r->[0]\n";
  }
}

print "got " . scalar( values( %H )) . " genes\n";



# read in (and display) expression matrix
$ta->loadFile($ARGV[1]);
my $a_ref2 = $ta->getArray();
my $r = shift @$a_ref2;
print join("\t", @$r); print "\n";
foreach my $r (@$a_ref2) {
  
  if (defined($H{ $r->[0] })) {
    $r->[0] = $H{ $r->[0] };
    print join("\t", @$r); print "\n";
  } else {
    #print "$r->[0] not there\n";
  }

}

