#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;
print join("\t", @$r) . "\n";

my %H = ();
foreach my $r (@$a_ref) {
  my $txt = join("\t", @$r);
  if (!defined($H{$txt})) {
    print "$txt\n";
    $H{$txt} = 1;
  }
}

