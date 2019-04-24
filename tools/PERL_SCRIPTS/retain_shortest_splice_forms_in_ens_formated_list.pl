BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %T_LEN = ();
my %T_RET = ();
foreach my $r (@$a_ref) {
  
  my $tlen = abs($r->[6] - $r->[5]);

  if (!defined($T_LEN{$r->[0]}) || ($tlen < $T_LEN{$r->[0]})) {
    $T_RET{$r->[0]} = $r;
    $T_LEN{$r->[0]} = $tlen;
  }
  
}


foreach my $r (values(%T_RET)) {
  print join("\t", @$r) . "\n";
}
