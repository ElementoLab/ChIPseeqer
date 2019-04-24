BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %CHR = ();
foreach my $r (@$a_ref) {
  
  push @{ $CHR{ $r->[ 1 ] } }, $r
  
}


foreach my $c (keys(%CHR)) {
  
  my $a_ref_c = $CHR{ $c };

  # order genes according to start/end
  
  my @a_tmp = sort { $a->[2] <=> $b->[2] } @$a_ref_c;

  for (my $i=0; $i<@a_tmp-1; $i++) {
    #print  . " -> " . ] . "\n";
    
    my $g1 = $a_tmp[$i]->[0];
    my $g2 = $a_tmp[$i+1]->[0];
    my $m = $a_tmp[$i+1]->[2] - $a_tmp[$i]->[3];
    
    print "$m\t$g1\t$g2\n";
    
  }
  
}

