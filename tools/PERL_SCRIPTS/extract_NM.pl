#!/usr/bin/perl


open IN, $ARGV[0] or die "Cannot open file";
my %H = ();
while (my $l = <IN>) {
  chomp $l;
  #my @a = split /\t/, $l, -1;
  
  while ($l =~ /(NM\_\d+)/g) {
    #print "$1\n";
    $H{$1} = 1;
  }


}
close IN;

print join("\n", keys(%H)) . "\n";
