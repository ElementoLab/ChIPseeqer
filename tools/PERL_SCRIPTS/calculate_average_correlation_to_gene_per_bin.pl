BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

if (@ARGV == 0) {
  die "args: bin_profile expression_matrix gene\n";
}


my $ta = Table->new;

$ta->loadFile($ARGV[1]);

#
# load expression matrix
#
my $a_ref_m = $ta->getArray();
my $h_ref_m = $ta->getIndexShifted();
my $pref    = $h_ref_m->{ $ARGV[2] }; die "Gene not found.\n" if !defined($pref);

my %COR     = ();
foreach my $r (@$a_ref_m) {
  
  my $n = shift @$r;
  next if ($n eq $ARGV[2]);

  my $p = Sets::pearson($pref, $r);
  
  $COR{$n} = $p;

}

#
# load bin profile
#
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
shift @$a_ref;

my %H = ();
foreach my $r (@$a_ref) {
  if (defined($COR{ $r->[0] })) {
    push @{ $H{ $r->[1] } }, $COR{ $r->[0] };
  } else {
    print STDERR "Attention, $r->[0] in expfile.\n";
  }
}

foreach my $g (sort(keys(%H))) {
  my $a = Sets::average($H{$g});
  my $s = Sets::stddev($H{$g});
  my $n = scalar( @{$H{$g}} );
  print sprintf("$g\t$n\t%4.3f\t%4.3f\n", $a, $s);

  #if (defined($ARGV[2])) {
   # Sets::writeSet($H{$g}, "$ARGV[2]/$g.txt");
  #}

}

