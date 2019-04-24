BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $h = shift @$a_ref;
print join("\t", @$h); print "\n";
shift @$h;

my %H = ();
foreach my $r (@$a_ref) {
  my $cg = shift @$r; 
  for (my $i=0; $i<@$r; $i++) {
    push @{ $H{ $cg }[ $i ] }, $r->[$i];
  }
}

foreach my $k (keys(%H)) {
  
  print "$k\t";
  my $r = $H{ $k };
  
  for (my $i=0; $i<@$r; $i++) {
    my $s = $H{$k}[$i];
    
    my %II = ();
    my $m  = @$s;
    die "$m cannot be 0\n" if ($m == 0);
    foreach my $t (@$s) {
      $II{ $t } ++;
    }

    if ($II{'P'} > $m/2) {
      print "\tP($II{P}/$m)";
    } else {
      print "\tA{$II{P}/$m)";
    }
    
  }	
  
  print "\n";
}

