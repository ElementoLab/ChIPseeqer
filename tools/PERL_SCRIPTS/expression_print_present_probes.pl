BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $cnt = 0;
foreach my $r (@$a_ref) {

  if ($cnt == 0) {
    print "$r->[0]\t$r->[1]\n";
  } else {
    
    if (($r->[2] eq "P") || ($r->[2] eq "M"))  {
      print "$r->[0]\t$r->[1]\n";
    }
    
  }

  

  $cnt ++;
}





