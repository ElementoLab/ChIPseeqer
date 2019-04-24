#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Motif;

use Getopt::Long;

if (@ARGV == 0) {
  print "Args: --motif=FILE --trimtype=[core,edges]\n";
}
my $motif    = undef;
my $trimtype = "edges";
my $te       = 1.25;
my $write    = 0;
my $outdir   = undef;

GetOptions("motif=s"    => \$motif,
           "trimtype=s" => \$trimtype,
	   "outdir=s"   => \$outdir,
	   "write=s"    => \$write,
	   "te=s"       => \$te);

my $mo = Motif->new;

$mo->readBulykWM($motif);

#$mo->print;
#$mo->motifEntropy();

if ($trimtype eq "core") {
  $mo->trimMotifToMostInformativeKmer();
} else { 
  $mo->shrinkMotifBasedOnEntropyV2($te);
}

if ($write == 0) {
  $mo->printBulyk;
} else {
  my $file = Sets::filename($motif);
  my @a = split /\_/, $file;
  my $co = $mo->getConsensus();
  #print "Consensus = $co\n";
  $a[0] = "TPBM";
  
  $a[1] = $co;
  shift @a;
  my $ff =  join("_", @a);
  if (defined($outdir)) {
    $ff = "$outdir/$ff";
  }
  if (-e $ff) {
    print "$ff already exist. Overwrite ?\n";
    <STDIN>;
  }
  open OUT, ">$ff" or die "Cannot open $ff\n";
  print OUT $mo->getBulykWMText();
  close OUT;	
  print "Created $ff\n";
}
