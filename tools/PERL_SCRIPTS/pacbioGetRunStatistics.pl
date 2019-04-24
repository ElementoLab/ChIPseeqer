#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $cntseq = 0;
my $avglen = 0;
my $maxlen = 0;
my $totlen = 0;
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    my $len = length($s);
    $avglen += $len;
    if ($len > $maxlen) {
      $maxlen = $len;
    }
    $cntseq++;
}
$totlen = $avglen;
$totlen = sprintf("%3.1f", $totlen/1000000);
$avglen = int(0.5+$avglen/$cntseq);
print "$ARGV[0]\t$cntseq\t$avglen\t$maxlen\t$totlen\n";
