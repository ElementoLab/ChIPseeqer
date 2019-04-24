#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Motif;

my $mo = Motif->new;

$mo->readBulykWM($ARGV[0]);

#$mo->print;
#$mo->motifEntropy();
$mo->shrinkMotifBasedOnEntropy($ARGV[1]);

if ($ARGV[2] eq "") {
  $mo->printBulyk;
} else {
  my $file = Sets::filename($ARGV[0]);
  my @a = split /\_/, $file;
  my $co = $mo->getConsensus();
  print "Consensus = $co\n";
  $a[1] = $co;
  my $ff = join("_", @a);
  if (-e $ff) {
    print "$ff already exist. Overwrite ?\n";
    <STDIN>;
  }
  open OUT, ">$ff" or die "Cannot open $ff\n";
  print OUT $mo->getBulykWMText();
  close OUT;	
  print "Created $ff\n";
}
