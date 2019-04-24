#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
if (@ARGV == 0) {
  die "Args: profile expfile\n";
}
my $ta = Table->new;

$ta->loadFile($ARGV[1]);

#
# load expression value
#
my $h_ref = $ta->getIndexKV(0,1);

#
# load 
#
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
shift @$a_ref;

my %H = ();
foreach my $r (@$a_ref) {
  if (defined($h_ref->{ $r->[0] })) {
    push @{ $H{ $r->[1] } }, $h_ref->{ $r->[0] };
  } else {
    #print STDERR "Attention, $r->[0] in expfile.\n";
  }
}

foreach my $g (sort(keys(%H))) {
  my $a = Sets::average($H{$g});
  #my $a = Sets::median($H{$g});
  my $s = Sets::stddev($H{$g});
  my $n = scalar( @{$H{$g}} );
  print sprintf("$g\t$n\t%4.3f\t%4.3f\n", $a, $s);

  if (defined($ARGV[2])) {
    Sets::writeSet($H{$g}, "$ARGV[2]/$g.txt");
  }

}

