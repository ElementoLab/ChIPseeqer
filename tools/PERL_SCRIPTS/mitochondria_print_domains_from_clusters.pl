BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

#
# read genes / domains
#
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndex(0);
my $a_ref_do  = $ta->getColumn(1);
$a_ref_do = Sets::removeDuplicates($a_ref_do);

#
# read genes 
#
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();
print "\t" . join("\t", @$a_ref_do) . "\n";
foreach my $r (@$a_ref) {
  
  #print "$r->[0]\n";
  
  my @a = split / /, $r->[0];
  my %CO = ();
  foreach my $o (@a) {
    next if ($o eq "");
    $CO{ $h_ref->{$o}->[1] }{ $h_ref->{$o}->[2] } ++;
  }
  
  my @b = split /\_/, $a[0];
  print "$b[0]";
  foreach my $k (@$a_ref_do) {
    my $n = scalar ( keys( % { $CO{ $k } } ) );
    print "\t$n";
  }
  print "\n";
  
}

