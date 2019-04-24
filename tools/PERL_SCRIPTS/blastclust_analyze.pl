BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->setDelim(" ");
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {

  my $n = @$r;

  my %H = ();
  my @a = ();
  foreach my $s (@$r) {
    my ($ss) = $s =~ /^(.+?)\_/;
    push @a, $ss;;
  }

  my $a_c = Sets::countSymbols(\@a);
  
  print "$n\t";
  foreach my $f (reverse(@$a_c)) {
    print "$f->[0]:$f->[1] " if ($f->[0] ne "");
  }

  print "\n";

  #print "$n\t" . join(" ", @a[0..10]) . "\n";


}

