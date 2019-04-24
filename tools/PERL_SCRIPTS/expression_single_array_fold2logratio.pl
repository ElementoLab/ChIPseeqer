BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $cnt = 0;
foreach my $r (@$a_ref) {
  if ($cnt == 0) {
    print join("\t", @$r); print "\n";
  } else {
    
    if ($r->[1] < 0) {
      $r->[1] = - 1 / $r->[1];
    } 

    print "$r->[0]\t" . sprintf("%5.4f", log($r->[1])/log(2)) . "\n";
    
  }
  
  $cnt ++;
}

