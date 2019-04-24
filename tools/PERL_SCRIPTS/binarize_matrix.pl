open IN, $ARGV[0];
my $l = <IN>;
print $l;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  my $n = shift @a; print "$n\t";
  
  foreach my $r (@a) {
    if ($r > 0) {
      $r = 1;
    }
  }
  
  print join("\t", @a); print "\n";
  
}

close IN;
