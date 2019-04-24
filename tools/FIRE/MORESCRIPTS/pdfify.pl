my $out = `find . -name \"*.eps\"`;

my @a = split /[\r\n]/, $out;

foreach my $f (@a) {
  print "$f .. ";
  my $fpdf = $f; $fpdf =~ s/\.eps$/\.pdf/;
  system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $f $fpdf");
  print "Done\n";
}
