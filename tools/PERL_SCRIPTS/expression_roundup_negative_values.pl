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
    
    my $g = shift @$r;

    foreach my $s (@$r) {
      $s = 1 if ($s < 0); 
    }
    
    print "$g\t" . join("\t", @$r); print "\n";
    
  }
  
  $cnt ++;
}

