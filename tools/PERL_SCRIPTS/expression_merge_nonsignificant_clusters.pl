#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;

my $a_ref = Sets::readSet($ARGV[0]);
my $i = 1;
my %H = ();
foreach my $r (@$a_ref) {
  $H{ $r } = $i++;
}

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();
my $r = shift @$a_ref;
print "$r->[0]\t$r->[1]\n";
foreach my $r (@$a_ref) {
  my $idx = 0;
  if (defined($H{$r->[1]})) {
    $idx = $H{$r->[1]};
  }	
  print "$r->[0]\t$idx\n";
}

