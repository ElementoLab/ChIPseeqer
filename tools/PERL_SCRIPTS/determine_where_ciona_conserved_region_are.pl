use lib qw(/home/elemento/PERL_MODULES);
use Sets;
use strict;

# read in GFF file
my @EXONS = ();
open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;
  if ($a[2] eq "exon") {
    my @a_tmp = (lc($a[0]), $a[3], $a[4], $a[8]);
    push @EXONS, \@a_tmp;
  }
}
close IN;


while (my $n = <STDIN>) {
  chomp $n;
  
  my @a = split / /, $n, -1;

  my $f = 0;
  foreach my $e (@EXONS) {
    if (($a[1] eq $e->[0]) && Sets::sequencesOverlap($e->[1], $e->[2], $a[2], $a[3])) {
      #print "$a[0] belongs to $e->[3]\n";
      $f = 1;
    } 
  }

  if ($f == 0) {
    print "$n\n";
  }
  
  
}
