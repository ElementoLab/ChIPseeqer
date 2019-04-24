open IN, $ARGV[0];
my $l = <IN>;
print "$l";
while (my $l = <IN>) {
  $l =~ s/\ +//g;
  print $l;

}
close IN;

