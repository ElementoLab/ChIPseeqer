BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $cnt = 0;
foreach my $r (@$a_ref) {
  if ($cnt == 0) {
    #print join("\t", @$r); print "\n";
    print "\thighest\n";
  }
  
  my $g = shift @$r;
  
  my $n = scalar(@$r);
  
  my $ii = undef;
  for (my $i=0; $i<$n; $i++) {
    if ($r->[$i] == 1) {
      $ii = $i;
    }
  }

  if (defined($ii)) {
    print "$g\t$ii\n";
  } 

  $cnt ++;
  
}

