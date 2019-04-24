#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Fasta;
use strict;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

print "GENE\tGC\tCpG\n";

while ( my $a_ref = $fa->nextSeq ) {
  my ($name, $seq) = @{$a_ref};
  
  $seq = uc($seq);
  
  my $l = length($seq);
  
  my @a = split //, $seq;
  
  my $gc  = 0;
  my $cpg = 0;

  my $effl = 0;  
  for (my $i=0; $i<@a; $i++) {
    if ($a[$i] ne 'N') {
      $effl ++;
    }
 
    if (($a[$i] eq 'C') || ($a[$i] eq 'G')) {
      $gc ++;
    }
    
    if (($a[$i] eq 'C') && ($a[$i+1] eq 'G')) {
      $cpg ++;
    }
    
  }
  
  $gc = $gc / $effl;
  $cpg = $cpg / ($effl-1);
  
  print "$name\t$gc\t$cpg\n";
    
}

