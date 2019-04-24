#BEGIN{ $home = `echo \$HOME`; chomp $home}
#use lib "$home/PERL_MODULES";

use lib qw(/home/elemento/PERL_MODULES);

use Sets;

my $ref = shift @ARGV;
my $a_ref = Sets::readSet($ref);

my %H = ();
my $c = 1;
foreach my $f (@ARGV) {
  my $a_ref_s = Sets::readSet($f);
  foreach my $g (@$a_ref_s) {
    $H{ $g } = $c;
  }
  $c ++;
}


foreach my $g (@$a_ref) {
  if (defined($H{ $g })) {
    print "$g\t$H{$g}\n";
  } else {
    print "$g\t0\n";
  }
  
}
