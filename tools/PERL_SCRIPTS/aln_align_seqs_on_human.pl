#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use ClustalW;
use Sets;

use strict;

my $cl = ClustalW->new;

#my $g = length("Orangutan           ");
#$cl->setNumCharName($g);

$cl->setFile($ARGV[0]);

my $a_ref_species = Sets::readSet("$ENV{HOME}/DATA/PROTEINS/ALN/species.txt");

my $a_ref_aln = $cl->getSeqsWithNames();

my %H = ();
foreach my $s (@$a_ref_aln) {
  my @a = split //, $s->[1];
  $H{$s->[0]} = \@a;
}

my $ah = $H{"Human"};

my $l  = @$ah;

my %NH = ();
my $el = 0;
for (my $i=0; $i<$l; $i++) {
  
  if ($ah->[$i] ne '-') {
    
    foreach my $s (@$a_ref_species) {
      if (defined($H{$s})) {
	push @{$NH{$s}}, $H{$s}->[$i];
      } else {
	push @{$NH{$s}}, '-';
      }
      
    }
    $el ++;
  }
} 

#foreach my $s (@$a_ref_species) {
#  print "$s\t$NH{$s}\n";
#}

for (my $i=0; $i<$el; $i++) {
  my $ii = $i+1;
  #print "$ii\t";
  foreach my $s (@$a_ref_species) {
    print $NH{$s}->[$i];
  }
  print "\n";
}
