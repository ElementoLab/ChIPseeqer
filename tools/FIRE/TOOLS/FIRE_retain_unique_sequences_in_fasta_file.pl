#!/usr/bin/perl


use lib "$ENV{FIREDIR}/SCRIPTS";

use Sets;

if ((!$ENV{FIREDIR}) || ($ENV{FIREDIR} eq '')) {
  die "Please ser the FIREDIR environment variable.\n";
}

use Fasta;

my %H = ();


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  if (!defined($H{$n})) {
    print ">$n\n$s\n\n";
    $H{$n} = 1;
  }
}

