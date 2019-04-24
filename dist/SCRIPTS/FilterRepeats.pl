open IN, $ARGV[0] or die "Cannot open file";

my $numrepeats = 0;

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  #my @b = split /\ /, $a[0];
  #if (($b[1] =~ /LINE/) || ($b[1] =~ /SINE/) || ($b[1] =~/DNA/)) {
    print "$a[0]-$a[1]-$a[2]-$numrepeats\t$a[3]\t$a[4]\t$a[5]\n";
  #}
  $numrepeats++;
}
close IN;

