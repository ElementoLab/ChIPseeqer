#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use Motif;
use strict;

my $scriptdir = "$ENV{FIREDIR}/SCRIPTS";


my $mo = Motif->new;

my $n = 10000;
$mo->readBulykWM($ARGV[0]);

my $txt = "";
for (my $i=0; $i<$n; $i++) {
  my $s = $mo->generateRandomSequence();
  $txt .= "$s\n";
}

my $tmpfile = Sets::getTempFile("/tmp/wm");

open OUT, ">$tmpfile" or die "Cannot opem $tmpfile\n";
print OUT "$txt";
close OUT;

my $w = int(0.5 + ($mo->getSize() * 5) / 9) * 1.5;

system("$scriptdir/weblogo/seqlogo -f $tmpfile -F EPS  -a -c -M -n -Y -w $w -h 3 > $ARGV[0].eps");

system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $ARGV[0].eps $ARGV[0].pdf");

system("open $ARGV[0].pdf");

unlink $tmpfile;






