#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;
print     Sets::jointab($r);

foreach my $r (@$a_ref) {
  my $n = shift @$r;
  
  print "$n";
  $s =~ s/\+//;
  foreach my $s (@$r) {
    if ($s > 0) {
      $s = -1 * $s;
    } 
    my $pvtext = Sets::round_pvalue_up(Sets::logp_to_p($s));
    if ( ($s > Sets::log10(0.01)) && ($ARGV[1] ne "") ) {
      print "\tNS";
    } else {
      print "\t$pvtext";
    }
  }
  
  print "\n";
  
}



