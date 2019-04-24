#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

# load refLink
my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndexKV(2,0);

# load index
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


my %H = ();
foreach my $r (@$a_ref) {

  if (defined($h_ref->{$r->[0]})) {
    $r->[0] = $h_ref->{$r->[0]};
  } else {
    print STDERR "Could not find match for $r->[0]\n";
  }
  my $n = shift @$r;
  push @{ $H{$n} }, @$r;

}

foreach my $k (keys(%H)) {
  my $r = Sets::removeDuplicates($H{$k});
  print "$k\t" . join("\t", @$r) . "\n";
}
