my $out = `find . -name \"*.eps\"`;

my @a = split /[\r\n]/, $out;

foreach my $f (@a) {
  print "$f .. ";
  my $fpdf = $f; $fpdf =~ s/\.eps$/\.png/;
  system("pstoimg -antialias $f");
  print "Done\n";
}
