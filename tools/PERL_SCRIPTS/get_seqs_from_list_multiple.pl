#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sequence;

use Sets;
use Table;
use strict;
#
#
# INPUT SP1, SP2, ORT, MIN, MAX
#

if (scalar(@ARGV) < 3) {
    die "Usage : getseqsfromlist.pl ORT DB1 DB2 DB3 .etc.\n";
}


my $list = shift @ARGV; 

foreach my $f (@ARGV) {
  unlink "$f.new" if (-e "$f.new");
}  

open INCH, $list;
my $cnt = 0;
while (my $l = <INCH>) {
  
  chomp $l;
  
  my $se = Sequence->new;
  
  my @SEQ = ();
  foreach my $r (@ARGV) {
    $se->setBlastDB($r);
    my $seq = $se->getSequenceFromBlastDB($l, 0, 0);
    
    if ($seq) {
      push @SEQ, $seq;
    }
  }
  
  if (scalar(@SEQ) == scalar(@ARGV)) {
    
    
    my $i = 0;
    foreach my $r (@ARGV) {
      open OUT, ">>$r.new";
      print OUT ">$l\n$SEQ[$i]\n\n";
      close OUT;
      $i++;
    }
  } 

  $cnt ++;

  print "done with $l\n";

  last if ($cnt == 5);
}

close INCH;
