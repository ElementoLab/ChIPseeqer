#BEGIN{ $home = `echo \$HOME`; chomp $home}
#use lib "$home/PERL_MODULES";

use lib qw(/home/elemento/PERL_MODULES);
use Table;
use Sets;

my $ta = Table->new;

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $nbE = $ARGV[1];
my $nbC = $ARGV[2];
my $cnt = 0;
foreach my $r (@$a_ref) {

  if ($cnt == 0) {
    print "\tRATIO\n";
  } else {
    
    my $n = shift @$r;
    
    my @v1 = ();
    my $p1 = 0;
    for ($i=0, $j=0; $i<$nbE; $i++, $j+=2) {
      if ($r->[$j+1] eq "Present") {
	$p1 ++;
      }
      push @v1, $r->[$j];
    }
    
    
    my @v2 = ();
    my $p2 = 0;
    for ($i=0, $j=2*$nbE; $i<$nbC; $i++, $j+=2) {
      if ($r->[$j+1] eq "Present") {
	$p2 ++;
      }
      push @v2, $r->[$j];
    }

    
    if (($p1 > 0.5 * $nbE) || ($p2 > 0.5 * $nbC)) {
      my $a1 = Sets::average( \@v1 );
      my $a2 = Sets::average( \@v2 );
      
      my $ra = log( $a1 / $a2 ) / log(2.0);
      $ra    = sprintf( "%4.3f", $ra);
      print "$n\t$ra\n";
    }
    
  }

  $cnt++;
  
}

