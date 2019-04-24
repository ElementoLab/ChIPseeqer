use lib qw(/home/elemento/PERL_MODULES);
use Sets;
use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref_zm = $ta->getIndexKV(0,1);


my $zt = 0;
my $mt = 0;
foreach my $v (values(%$h_ref_zm)) {
  if ($v eq "zygotic") {
    $zt ++;
  }
  if ($v eq "maternal") {
    $mt ++;
  }
  
}

my $h_ref = Sets::getDictionnary($ARGV[0]);

# read in zygotic/maternal



foreach my $k (keys(%$h_ref)) {
  
  my $z = 0;
  my $m = 0;
  foreach my $g (@{ $h_ref->{ $k } }) {
    
    if (defined($h_ref_zm->{ $g })) {

      if ($h_ref_zm->{ $g } eq "zygotic") {
	$z ++;
      }
      if ($h_ref_zm->{ $g } eq "maternal") {
	$m ++;
      }
      
      
    }
    
  }

  my $ra = $mt/( $zt+$mt);

  my $ra0 = $m / ($z + $m);
  my $exp = ($z+$m);
  print "$k\t$m\t$exp\t$ra\n";
  
}
