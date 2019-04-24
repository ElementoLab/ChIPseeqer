BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %H1 = ();
my %H2 = ();
foreach my $r (@$a_ref) {
  next if ($r->[1] eq ""); 
  $H1{ $r->[0] } ++;
  $H2{ $r->[1] } ++;
  
}


foreach my $r (@$a_ref) {
  next if ($r->[1] eq ""); 
 
  if (($H1{ $r->[0] } == 1) && ($H2{ $r->[1] } == 1)) {
    print "$r->[0]\t$r->[1]\n";
  }
  
}



