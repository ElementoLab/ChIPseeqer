BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;



#
# load windows
#

use Table;
$ta->loadFile($ARGV[1]);
my $a_ref_w = $ta->getArray();



#
# load genes
#
$ta->loadFile($ARGV[0]);
my $a_ref_g = $ta->getArray();

foreach my $g (@$a_ref_g) {
  
  my $yes = 0;
  foreach my $w (@$a_ref_w) {
  
    next if ($g->[1] ne $w->[0]);  # same chr only

    # calculate distance with st and en
    
    my $ds = abs( $w->[1] - $g->[2] );
    my $de = abs( $w->[1] - $g->[3] );
    
    my $d  = Sets::min($ds, $de);
    if ($d < 5000) {
      $yes = $w->[2]; 
    }
  }
  
  if ($yes > 0) {
    print "$g->[0]\t$yes\n";
  }
  
}

