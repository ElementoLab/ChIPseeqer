use Getopt::Long;

my $peakfile = undef;
my $norm     = "rpkm";

if (@ARGV == 0) {
  die "--peakfile=FILE --norm=rpkm|numnt\n";
}


GetOptions("peakfile=s" => \$peakfile,
           "norm=s"     => \$norm);

#my $peakfile = $ARGV[0];

my $normtxt = "";
if ($norm eq "rpkm") {
  $normtxt = " -rpkmnorm 1 ";
} elsif ($norm eq "numnt") {
  $normtxt = " -numntnorm 1 ";
}


# EXTRACT 2kb AROUND CENTERS
my $todo = "perl $ENV{CHIPSEEQERDIR}/SCRIPTS/extract_regions_around_peak_summits.pl --peakfile=$peakfile --w=2000 --force=1 ";
system($todo);

# NUCLEOSOMES CB
$todo = "$ENV{CHIPSEEQERDIR}/ChIPseeqerGetReadDensityProfiles.bin -intervals $peakfile.centered2000 -chipdir $ENV{CHIPSEEQERDIR}/DATA/SOLEXA/Nucleosomes-NB -format bed -fraglen 0 $normtxt -uniquereads 0 -outfile $peakfile.centered2000.NUCLNB";
system($todo);

# NUCLEOSOMES LY1
$todo = "$ENV{CHIPSEEQERDIR}/ChIPseeqerGetReadDensityProfiles.bin -intervals $peakfile.centered2000 -chipdir $ENV{CHIPSEEQERDIR}/DATA/SOLEXA/Nucleosomes-LY1 -format bed -fraglen 0 $normtxt -uniquereads 0 -outfile $peakfile.centered2000.NUCLLY1";
system($todo);

# COMBINE
$todo = "expression_concatenate_matrices.pl $peakfile.centered2000.NUCLNB $peakfile.centered2000.NUCLLY1 > $peakfile.centered2000.NUCLNB.NUCLLY1";
system($todo);

# DRAW COMBINED AVERAFE
$todo = "R --slave --args $peakfile.centered2000.NUCLNB.NUCLLY1 < $ENV{CHIPSEEQERDIR}/SCRIPTS/makeAvgConsProfile.R";
system($todo);


# MAIL
$todo = "uuencode $peakfile.centered2000.NUCLNB.NUCLLY1 $peakfile.centered2000.NUCLNB.NUCLLY1 | mail ole2001\@med.cornell.edu";
system($todo);
