BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

my $m = $ARGV[0];

my ($gl) = $m =~ /^(\.*)/;
my ($gr) = $m =~ /(\.*)$/;

$m =~ s/^(\.*)//;
$m =~ s/(\.*)$//;


# print "$gl - $gr\n";

my $a_ref  = Sets::get_array_from_re($m);

my $a_shu  = Sets::shuffle_array($a_ref);

print $gl;

foreach my $r (@$a_shu) {
  print "$r";
}

print $gr;

print "\n";

