BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  my @A = ();

  for (my $i=0; $i<=$ARGV[2]; $i++) {
    push @A, $r->[$i];
  }
  push @A, $ARGV[1];

  for (my $i=$ARGV[2]+1; $i<@$r; $i++) {
    push @A, $r->[$i];
  }
  
  print join("\t", @A); print "\n";


}

