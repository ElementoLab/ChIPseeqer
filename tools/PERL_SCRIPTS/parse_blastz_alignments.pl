use lib qw(/home/elemento/PERL_MODULES);
use Sets;
use strict;

while (1) {

  my $n = <STDIN>;
  if (!$n) {
    exit;
  }
  my $s1 = <STDIN>; chomp $s1;
  my $s2 = <STDIN>; chomp $s2;
  my $l  = <STDIN>;

  next if (length($s1) < 100);

  my $i = Sets::getSequencesIdentity($s1, $s2);
  
  next if ($i < (defined($ARGV[0])?$ARGV[0]:0.7));

  print $n;
  
}
