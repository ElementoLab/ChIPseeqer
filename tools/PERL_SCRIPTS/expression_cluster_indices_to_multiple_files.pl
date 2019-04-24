use strict;

open IN, $ARGV[0] or die "cannot open $ARGV[0]\n";
my $cnt = 0;
my $h = <IN>;
my %H = ();
my %K = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;
  $H{ $a[0] } =  $a[1];
  $K{ $a[1] } =  1;
  
}
close IN;

foreach my $k (keys(%K)) {
  open OUT, ">$ARGV[0].c$k";
  print OUT $h;
  foreach my $g (keys(%H)) {
    if ($H{$g} == $k) {
      print OUT "$g\t1\n";
    } else {
      print OUT "$g\t0\n";
    }
  }
  close OUT;

}
