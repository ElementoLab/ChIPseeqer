use strict;

#
# load NM to gene name
#
my $refgene = "/Users/olivier/PROGRAMS/CHIP-SEQ/DATA/refGene.txt.7June2009";
my $h_ref   = {};
open IN, $refgene or die "Cannot open $refgene\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $h_ref->{$a[1]} = $a[12]
}
close IN;



open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";

my %IDX = ();

while (my $l = <IN>) {
  chomp $l;
  next if ($l =~ /\#/);
  next if ($l eq "");

  if ($l =~ /^a/) {
    my ($sc) = $l =~ /score\=(\d+)/;

    my $l1   = <IN>;
    chomp $l1;
    my @a1 = split /\ +/, $l1, -1;

    my $l2   = <IN>;
    chomp $l2;
    my @a2 = split /\ +/, $l2, -1;
    
    next if ($a2[4] eq '-');

    my @a_tmp = ($a1[1], $sc, $a1[2], $a1[3]);

    push @{ $IDX{ $a2[1] } }, \@a_tmp;

    #print "a=$l\ns=$l1\ns=$l2\n\n";
  }
}
close IN;

foreach my $i (keys(%IDX)) {

  print "$i";

  my @a_sorted = sort { $b->[1] <=> $a->[1] } @{$IDX{$i}};


  # move while genename and scores are same
  my $n = @a_sorted;
  my $g = $h_ref->{$a_sorted[0]->[0]};
  my $s = $a_sorted[0]->[1];

  my $j = 1;
  while (($j < $n) && (($g ne "") && ($g eq $h_ref->{$a_sorted[$j]->[0]})) && ($s == $a_sorted[$j]->[1])) {
    $j++;
  }

  if ($j == $n) {
    print "\tKEEP(1)";
  } else {
    # if previous score greater than current, keep
    if ($a_sorted[$j-1]->[1] > $a_sorted[$j]->[1]) {
      print "\tKEEP(2)";
    } else {
      print "\tSKIP";
    }
  }

  print "\t$j";

  for (my $k=0; $k<$j; $k++) {
    my $r = $a_sorted[$k];
    print "\t$r->[0]($h_ref->{$r->[0]})\t$r->[1]\t$r->[2]\t$r->[3]";
  }

  print "\n";

}
