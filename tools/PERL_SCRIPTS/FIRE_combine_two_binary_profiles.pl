#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
#shift @$a_ref;

my %H = ();
foreach my $r (@$a_ref) {
  my $txt = Sets::tabjoin($r);
  $H{ $txt } = 1;
  print $txt;
}

$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();
shift @$a_ref;

foreach my $r (@$a_ref) {
  my $txt = Sets::tabjoin($r);
  if (!defined($H{$txt})) {
    print $txt;
  }
}

