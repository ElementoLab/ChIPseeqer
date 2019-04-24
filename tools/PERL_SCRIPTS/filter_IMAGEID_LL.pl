BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  if (($r->[1] ne "") && ($r->[1] !~ /\;/) && ($r->[1] !~ /cluster/) && ($r->[1] !~ /Data/)) {
    print "$r->[0]\t$r->[1]\n";
  } 
}

