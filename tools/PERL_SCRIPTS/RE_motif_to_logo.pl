#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use Motif;

my $scriptdir = "$ENV{FIREDIR}/SCRIPTS";

my $txt = Sets::myre2scanace($ARGV[0]);
my $m   = Sets::re_num_positions($ARGV[0]);

my $tmpfile = Sets::getTempFile("/tmp/wm");

open OUT, ">$tmpfile" or die "Cannot opem $tmpfile\n";
print OUT "$txt";
close OUT;

my $w = int(0.5 + ($m * 5) / 9);

system("$scriptdir/weblogo/seqlogo -f $tmpfile -F EPS  -a -c -M -n -Y -w $w -h 3 > $ARGV[0].eps");

system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $ARGV[0].eps $ARGV[0].pdf");

unlink $tmpfile;






