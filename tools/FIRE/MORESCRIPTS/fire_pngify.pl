
my $out = `find $ARGV[0]\_FIRE/DNA_RNA/ -name \"*.eps\"`;

my @a = split /[\r\n]/, $out;

foreach my $f (@a) {
  my $fpng = $f; $fpng =~ s/\.eps$/\.png/;

  next if ($f =~ /\d.eps$/);
  next if (-e $fpng);

  print "$f .. ";
  system("pstoimg -antialias $f");

  my $fpdf = $f; $fpdf =~ s/\.eps$/\.pdf/;
  #system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $f $fpdf");

  print "Done\n";
}
