#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use Motif;
use strict;

my $scriptdir = "$ENV{FIREDIR}/SCRIPTS";

my $txt = Sets::myre2wm_aa($ARGV[0]);
my $m   = Sets::re_num_positions($ARGV[0]);

my $tmpfile = Sets::getTempFile("/tmp/wm");

open OUT, ">$tmpfile" or die "Cannot opem $tmpfile\n";
print OUT "$txt";
close OUT;

my $w = 3 * int(0.5 + ($m * 5) / 9);

print "Size = $w ($m pos)\n";

#system("$scriptdir/weblogo/seqlogo -f $tmpfile -F EPS -k 0 -S -a -c -M -n -Y -w 5 -h 3 > $ARGV[0].eps");

system("$scriptdir/weblogo/seqlogo -f $tmpfile -F EPS -k 0 -S -a -c -M -n -Y -w $w -h 3 > $ARGV[0].eps");

system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $ARGV[0].eps $ARGV[0].pdf");

unlink $tmpfile;






