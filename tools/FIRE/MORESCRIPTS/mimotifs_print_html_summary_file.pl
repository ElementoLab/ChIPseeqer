BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

my $a_ref = Sets::getFiles($ARGV[0]);

print "<table cellspacing=0>\n";
foreach my $r (@$a_ref) {

  my $f1 = undef;
  $f1 = "<a href=\"$r.$ARGV[1].pdf\">upstream</a>" if (-e "$r.$ARGV[1].pdf");
  my $f2 = undef;
  $f2 = "<a href=\"$r.$ARGV[2].pdf\">3'UTRs</a>" if (-e "$r.$ARGV[2].pdf");

  $r =~ s/\_/ /g;
  $r =~ s/\.txt\.nodups//;
  print "<tr><td><tt>$r</tt></td><td><tt>$f1&nbsp;</tt></td><td><tt>&nbsp;$f2</tt></td></tr>\n";
}
print "</table>\n";
