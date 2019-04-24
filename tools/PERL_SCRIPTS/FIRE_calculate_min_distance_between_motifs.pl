#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use strict;
my $ta = Table->new;
#srand(1234);

$ta->loadFile("/Users/olivier/PEOPLE/WEIMIN/DESIGN/OLIVIER/promoter_regions_lengths.txt");
my $h_ref_len = $ta->getIndexKV(0,1);



$ta->loadFile("$ARGV[0]");
my $h_ref_exp = $ta->getIndexKV(0,1);


$ta->loadFile("$ARGV[0].profiles");
my $a_ref = $ta->getArray();

my %H = ();
my %G = ();
foreach my $r (@$a_ref) {
  
  push @{ $H{ $r->[0] }{ $r->[1] } }, $r->[2];
  $G{$r->[1]} = 1;
}


my @motifs = keys(%H);
my @genes  = keys(%G);

#print scalar(@motifs) . " motifs\n";
#print scalar(@genes) . " genes\n";

for (my $i=0; $i<@motifs-1; $i++) {
  
  for (my $j=$i+1; $j<@motifs; $j++) {
    
    foreach my $g (@genes) {

      next if (!defined($H{$motifs[$i]}{$g}) || !defined($H{$motifs[$j]}{$g}));

      next if ($h_ref_exp->{$g} != 1);
      
      my @pos_gene1 = @{ $H{$motifs[$i]}{$g} };
      my @pos_gene2 = @{ $H{$motifs[$j]}{$g} };

      next if ((@pos_gene1 == 0) || (@pos_gene2 == 0));

      my $min_pos = 100000000;

      for (my $k=0; $k<@pos_gene1; $k++) {
	
	for (my $l=0; $l<@pos_gene2; $l++) {

	  my $p1 = $pos_gene1[$k];
	  my $p2 = $pos_gene2[$k];

	  if ($ARGV[1] ne "") {
	    my $p1 = int( 0.5 + rand() * $h_ref_len->{$g} ); 
	    my $p2 = int( 0.5 + rand() * $h_ref_len->{$g} );
	  }

	  my $d = abs($p1 - $p2);
	  
	  if ($d < $min_pos) {
	    $min_pos = $d;
	  }
	  
	}	
	

	
      }
      print "$g\t$min_pos\n";
      
    }
   
  }
  
}
