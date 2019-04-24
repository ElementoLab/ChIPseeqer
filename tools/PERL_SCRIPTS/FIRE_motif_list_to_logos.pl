#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $scriptdir = "$ENV{FIREDIR}/SCRIPTS";

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  
  print "processing $r->[0]\n";

  my $txt = Sets::myre2wm($r->[0]);
  
  my $tmpfile = Sets::getTempFile("/tmp/wm");
  
  open OUT, ">$tmpfile" or die "Cannot opem $tmpfile\n";
  print OUT "$txt\n";
  close OUT;

  system("$scriptdir/weblogo/seqlogo -f $tmpfile -F EPS  -a -c -M -n -Y -w 5 -h 3 > M$r->[0].eps");

  system("ps2pdf -dEPSCrop -dAutoRotatePages=/None M$r->[0].eps M$r->[0].pdf");

  unlink $tmpfile;

}




