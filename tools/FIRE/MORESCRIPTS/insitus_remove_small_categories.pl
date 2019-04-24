BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;

my $ta    = Table->new; 
my $a_ref = Sets::getFiles($ARGV[0]);

mkdir "SMALL" if (! -e 'SMALL');

foreach my $f (@$a_ref) {

  print "$f\n";
  $ta->loadFile($f);
  my $a_ref_t = $ta->getArray();
  shift @$a_ref_t;
  my $cnt = 0;
  foreach my $r (@$a_ref_t) {
    $cnt ++ if ($r->[1] == 1);
  }	
  
  print "got $cnt 1s.\n";
  
  if ($cnt < 100) {
    system("mv $f SMALL/$f");
  }

}
