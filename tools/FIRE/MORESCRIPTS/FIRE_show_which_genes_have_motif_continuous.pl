BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Table;
use Sets;

# load profiles
my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref_pro = $ta->getArray();

foreach my $r (@$a_ref_pro) {
  if ($r->[0] eq $ARGV[2]) {
    push @{ $H{ $r->[1] } }, $r->[2];
  }
}


# load expfile
$ta->loadFile($ARGV[0]);
my $a_ref_exp = $ta->getArray();

my %I = ();
foreach my $r (@$a_ref_exp) {
  print join("\t", @$r) . " ";
  if (defined($H{ $r->[0] })) {
    print "*";
    $I{$r->[0]} = 1;
  }
  print "\n";
}
  
# go thru profiles again
foreach my $r (@$a_ref_pro) {
  if (($r->[0] eq $ARGV[2]) && defined($I{$r->[1]})) {
    print join("\t", @$r) . "\n";
  }
}

