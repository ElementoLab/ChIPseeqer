BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $i = 0;
foreach my $r (@$a_ref) {
  shift @$r; 
  my $n = shift @$r;
  my $txt = "$n";
  my $P   = 0;
  for (my $i=0; $i<@$r; $i+=2) {
    $txt .=  "\t$r->[$i]";
    if ($r->[$i+1] eq 'P') {
      $P ++;
    }
  }
  
  if (($i == 0) || ($P >= 10)) {
    print "$txt\n";
  } 

  $i++;
}

