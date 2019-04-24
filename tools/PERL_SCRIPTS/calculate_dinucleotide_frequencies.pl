#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Fasta;
use strict;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $ltot = 0;
my %COUNT = ();
while ( my $a_ref = $fa->nextSeq ) {
    my ($name, $seq) = @{$a_ref};
 
    my $l = length($seq);
    $seq  = uc($seq);
    #$ltot += $l;
    my @a = split //, $seq;
    
    for (my $i=1; $i<@a; $i++) {
      my $nt1 = $a[$i-1];
      my $nt2 = $a[$i  ];
      if (($nt1 =~ /[ACGT]/) && ($nt2 =~ /[ACGT]/)) {
	$COUNT{"$nt1$nt2"} ++;
	$ltot ++;
      }
    }
}

foreach my $k (sort(keys(%COUNT))) {
  print "$k\t" . ($COUNT{$k}/$ltot) . "\n";
}



