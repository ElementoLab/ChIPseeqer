BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %CLU = ();
foreach my $r (@$a_ref) {

  push @{ $CLU{ $r->[1] } }, $r->[0];

}


foreach my $k (keys(%CLU)) {
  
  print "$k\n";
  
  open OUT, ">c$k.txt" or die "cannot open c$k.txt\n";
  
  print OUT join("\n", @{ $CLU{ $k } }) . "\n";
  
  close OUT;
  
}

