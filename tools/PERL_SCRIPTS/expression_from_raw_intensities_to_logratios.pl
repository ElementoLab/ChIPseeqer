BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $f = shift @ARGV;

my $ta = Table->new;
$ta->loadFile($f);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
my @b = ();
push @b, $r->[0];
for (my $i=0; $i<@ARGV; $i++) {    
  my @a = split /\//, $ARGV[$i];    
  my $s = $r->[ $a[0] ] . "/" . $r->[ $a[1] ];
  push @b, $s;
}
print join("\t", @b) . "\n";


foreach my $r (@$a_ref) {
  
  my @b = ();
  push @b, $r->[0];
  for (my $i=0; $i<@ARGV; $i++) {    
    my @a = split /\//, $ARGV[$i];    
    my $s = $r->[ $a[0] ] / $r->[ $a[1] ];
    $s    = sprintf("%4.3f", log($s)/log(2.0));
    push @b, $s;
  }
  print join("\t", @b) . "\n";
  
}

