BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  if ($r->[2] eq 'protein') {
    my ($n) = $r->[8] =~ /Alias\=(.+?)\;/;
    my $s = ($r->[6] eq '+'?1:-1); 
    print "$n\t$r->[0]\t$r->[3]\t$r->[4]\t$s\t$r->[3]\t$r->[4]\n";
  }	
}

