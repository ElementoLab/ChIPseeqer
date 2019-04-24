#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use ClustalW;
use strict;

my $cl = ClustalW->new;

my $g = 26;
if ($ARGV[1] ne "") {
  $g = $ARGV[1];
}

#$cl->setNumCharName($g);

$cl->setFile($ARGV[0]);

my $a_ref_aln = $cl->getSeqsWithNames();


my $n = @$a_ref_aln;


my @alns = ();
my @names = ();
my $i = 0;

for (my $i=0; $i<$n-1; $i++) {
  
  for (my $j=0; $j<$n-1; $j++) {
    
    my $d = Sets::getDistanceBetweenSequences($a_ref_aln->[$i]->[1], $a_ref_aln->[$j]->[1]); 
    print "$a_ref_aln->[$i]->[0]\t$a_ref_aln->[$i]->[0]\t$d\n";
    
  }

}

