BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  print "$r->[0]\t$r->[1]\n";

  my $re = Sets::consensus2re($r->[1]);
  my $a  = Sets::get_array_from_re($re);
  my $n  = @$a;
  open IN, ">ACE/$r->[0].txt";
  print IN "Motif 1\n";
  print IN Sets::myre2wm($re);
  print IN '*' x $n;
  close IN;
  
}

