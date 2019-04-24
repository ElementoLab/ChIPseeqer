#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();
shift @$a_ref;

my %H = ();
foreach my $r (@$a_ref) {
  #my $txt = Sets::tabjoin($r);
  $H{ $r->[0] } = $r->[1];
  #print $txt;
}

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
shift @$a_ref;

print "i\ti\n";
foreach my $r (@$a_ref) {
  #my $txt = Sets::tabjoin($r);
  #if (!defined($H{$r->})) {
  #  print $txt;
  #}
  
  if (($r->[1] == 1) || ($H{$r->[0]} == 1)) {
    print "$r->[0]\t1\n";
  } else {
    print "$r->[0]\t0\n";
  }

}

