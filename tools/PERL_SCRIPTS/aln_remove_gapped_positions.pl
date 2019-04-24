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



my @alns = ();
my @names = ();
my $i = 0;
foreach my $s (@$a_ref_aln) {  
  my @a = split //, $s->[1];
  $alns[$i] = \@a;
  push @names, $s->[0];
  $i++;
}

my $l = @{$alns[0]};
my $n = @alns;

my @newalns = ();
for (my $i=0; $i<$l; $i++) {
  my $gapped = 0;
  for (my $j=0; $j<$n; $j++) {
    #print "$alns[$j][$i]\n";
    if ($alns[$j][$i] eq '-') {
      $gapped = 1;
      last;
    }
  }  
  #print "gapped = $gapped\n";
  if ($gapped == 0) {
    for (my $j=0; $j<$n; $j++) {
      $newalns[$j] .=  $alns[$j][$i];
    }
  }
}


my $cl2 = ClustalW->new;
$cl2->setNumCharName($g);
$cl2->setSequences(\@names, \@newalns);
print $cl2->getClustalWformat;
