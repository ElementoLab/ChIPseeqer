BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
my $ta = Table->new;

if (@ARGV == 0) {

}

#
# /home/elemento/PEOPLE/WEIMIN/METHYLATION/PLATFORM/prom_coordinates.txt /home/elemento/PEOPLE/WEIMIN/DESIGN/OLIVIER/BCL6_bound_peak_regions.txt
#

#
# load transcripts
#
$ta->loadFile($ARGV[1]);
my $a_ref_transcripts = $ta->getArray();


# load regions
$ta->loadFile($ARGV[0]);


my $a_ref_regions = $ta->getArray();

foreach my $m (@$a_ref_regions) {

  my $ov = 0;
  foreach my $b (@$a_ref_transcripts) {
    
    next if ($m->[1] ne $b->[1]);

    my $tss = undef;
    if ($b->[4] == 1) {
      $tss = $b->[5];
    } else {
      $tss = $b->[6];
    }
    
    my $dmin = Sets::min( abs($m->[2] - $tss), abs($m->[3] - $tss) );

    if ($dmin < 500) {
      print "$m->[0]\t";
      #print "$m->[2]\t$m->[3]\t";
      print "$b->[0]";
      #print "\t";
      #print "$b->[5]\t$b->[6]\t";
      print "\n";
    }

  }

}


