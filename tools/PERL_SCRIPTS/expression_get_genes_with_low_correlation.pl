BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $h_ref = $ta->getIndexShifted();

my $p = $h_ref->{$ARGV[1]};

shift @$a_ref;

my %H = ();
foreach my $r (@$a_ref) {

  my @p1 = @{$r}[1..3]; 
  my @p2 = @{$r}[4..6];
  $H{ $r->[0] } = Sets::pearson(\@p1, \@p2);
 
}

my $o_ref = Sets::hash_order(\%H);
foreach my $k (reverse(@$o_ref)) {
  if ( $H{$k} < 0.75 ) {
    print "$k\t" . sprintf("%4.3f", $H{$k}) . "\n";
  }	
}

