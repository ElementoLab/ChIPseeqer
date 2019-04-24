BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref; print join("\t", @$r). "\n";

my @A = ();
my $i = 0;
foreach my $r (@$a_ref) {
  my @o = (@$r, $i);
  push @A, \@o;
  $i ++;
}


@A = sort { $a->[1] <=> $b->[1] } @A;

my $i = 0;
foreach my $r (@A) {
  print "$r->[0]\t$i\n";
  $i ++;
}

