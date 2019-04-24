foreach my $f (@ARGV) {
  if (-e $f) {
    system("cat $f");
  }
}
