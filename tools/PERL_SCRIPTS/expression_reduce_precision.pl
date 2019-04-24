BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;

print join("\t", @$r) . "\n";

foreach my $r (@$a_ref) {
  for (my $i=1; $i<@$r; $i++) {
    
    my $prec = "3.2";
    if ($ENV{PREC} ne "") {
      my $p1 = $ENV{$PREC}+1;
      $prec = "$p1.$ENV{PREC}";
    }
    $r->[$i] = sprintf("%$prec" . "f", $r->[$i]);
  }
  print join("\t", @$r) . "\n";
}

