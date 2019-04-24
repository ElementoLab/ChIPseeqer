#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use ClustalW;
use Sets;

use strict;

my $cl = ClustalW->new;

$cl->setFile($ARGV[0]);

my $a_ref_aln = $cl->getSeqsWithNames();

my $s1 = $a_ref_aln->[0]->[1];
my $s2 = $a_ref_aln->[1]->[1];

print "$s1\n";
print "$s2\n";

my @a1 = split //, $s1;
my @a2 = split //, $s2;

while ($a1[0] eq '-') {
  shift @a1;
  shift @a2;
}

#print join("", @a1) . "\n";

my $i1 = 0;
my $i2 = 0;
for (my $j=0; $j<@a1; $j++) {
  
  print "$a1[$j]$a2[$j]\n";

  if ($a1[$j] ne '-') {
    $i1++;
  }
  if ($a2[$j] ne '-') {
    $i2++;
  }

  if ($i1 == $ARGV[1]) {
    print "$i2 ($a2[$j])\n";
    last;
  }


}


