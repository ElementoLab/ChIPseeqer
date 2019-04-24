BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

my %CHR = ();
foreach my $r (@$a_ref) {
  push @{ $CHR{ $r->[0] } }, $r if (defined($r->[1]) && ($r->[1] ne ""));
}


open IN, $ARGV[0];

open OUT1, ">blastz_species1.fa" or die "Cannot open file 1\n";
open OUT2, ">blastz_species2.fa" or die "Cannot open file 2\n";;

my @a = ();
my $cnt = 0;
while (my $l = <IN>) {
  chomp $l;
  next if ($l =~ /^\#/);
  next if ($l eq "");

  if ($l =~ /^a\ /) {
    
    my ($sc) = $l =~ /score\=(.+)$/;


    next if ($sc < 2000); #want high scoring elements 

    my $l1 = <IN>;
    my @a1 = split /\ +/, $l1;
    my $l2 = <IN>;
    my @a2 = split /\ +/, $l2;
    
    next if ($a1[3] < 200);  # want elements of size >= 200

    my ($species, $chr) = split /\./, $a1[1];



    my $s  = $a1[2];
    my $e  = $s + $a1[3];
    
    my $ov = 0;
    
    foreach my $g (@{ $CHR{ $chr } }) {      
      if (Sets::sequencesOverlap($g->[1], $g->[2], $s, $e)) {
	$ov = 1; last;	
      }
    }
    
    if ($ov == 0) {

      print "$cnt.$a1[1].$s.$e\n";
      
      $a1[6] =~ s/\-//g;
      $a2[6] =~ s/\-//g;
      
      print OUT1 ">$cnt.$a1[1].$s.$e\n$a1[6]\n\n";
      print OUT2 ">$cnt.$a2[1]\n$a2[6]\n\n";
      
      $cnt ++;
    }
    
  }

  
  
}

close IN;
close OUT1;
close OUT2;
