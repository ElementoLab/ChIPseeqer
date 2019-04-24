#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;

my $a_ref_f = Sets::getFiles($ARGV[0]);

my $a_ref_AP2 = Sets::readSet($ARGV[1]);
foreach my $a (@$a_ref_AP2) {
  $a = Sets::Pf_motifname2genename($a);

}
foreach my $f (@$a_ref_f) {
  my $ff = $f;
  $ff =~ s/^.+varnorm\.//;
  my $reg = Sets::Pf_motifname2genename($ff);

  open IN, $f;
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    print "$reg\tpd\t$a[0]\n" if (($a[1] > 0) && Sets::in_array($a[0], @$a_ref_AP2));
  }
  close IN;



}
