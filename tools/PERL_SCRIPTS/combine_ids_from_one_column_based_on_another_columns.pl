BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %H = ();
foreach my $r (@$a_ref) {  
  next if ($r->[$ARGV[1]] eq "");
  next if ($r->[$ARGV[2]] eq "");
  push @{ $H{ $r->[$ARGV[1]] } }, $r->[$ARGV[2]];

}


foreach my $k (keys(%H)) {
  print "$k\t" . join("/", @{$H{$k}}) . "\n";
}



