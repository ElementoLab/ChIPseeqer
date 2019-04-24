#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;



use Table;

my $ta = Table->new;

# load first file
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

# load second file
$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndexKV(0,1);

shift @$a_ref;
print "GENE\tEXP\n";
foreach my $r (@$a_ref) {

  next if (!defined($h_ref->{$r->[0]}));

  my $i = "?";
  if (($r->[1] eq "BIVALENT") && ($h_ref->{ $r->[0] } eq "BIVALENT")) {
    $i = "H1:BIVALENT -> hVPrEC:BIVALENT";
  } elsif (($r->[1] eq "BIVALENT") && ($h_ref->{ $r->[0] } eq "H3K4me3")) {
    $i = "H1:BIVALENT -> hVPrEC:H3K4me3 (RESOLVED)";
  } elsif (($r->[1] eq "BIVALENT") && ($h_ref->{ $r->[0] } eq "H3K27me3")) {
    $i = "H1:BIVALENT -> hVPrEC:H3K27me3 (RESOLVED)";
  } elsif (($r->[1] eq "BIVALENT") && ($h_ref->{ $r->[0] } eq "NONE")) {
    $i = "H1:BIVALENT -> hVPrEC:NONE";
  } else {
    $i = "H1:$r->[1] -> hVPrEC:$h_ref->{$r->[0]}";
  }

  print "$r->[0]\t$i\n";

}
