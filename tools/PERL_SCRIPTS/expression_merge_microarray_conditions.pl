BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;

my %keys   = ();
my @tables = ();
foreach my $f (@ARGV) {
  my $ta = Table->new;
  $ta->loadFile($f);
  my $h_ref = $ta->getIndexKV(0,1);

  push @tables, $h_ref;

  my @tk = keys(%$h_ref);

  foreach my $r (@tk) {
    $keys{ $r} = 1;
  }

}

print "ID_REF";
foreach my $f (@ARGV) {
  print "\t$f";
}
print "\n";

foreach my $k (sort(keys(%keys))) {

  next if ($k eq "ID_REF");
  print $k;


  
  foreach my $t (@tables) {
    print "\t" . $t->{ $k }; 
  }
  print "\n";

}
