BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

for (my $i=0; $i<@$a_ref; $i+=2) {
  
  if (($a_ref->[$i]->[3] == 0) && ($a_ref->[$i+1]->[3] == 1)) {
    print join("\t", @{ $a_ref->[$i  ] } ) . "\n";
    print join("\t", @{ $a_ref->[$i+1] } ) . "\n";
  }
  
}

