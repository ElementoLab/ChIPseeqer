#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;
use PBS;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args --files=FILE --chipdir=DIR --format=STR\n";
}

my $chipdir = "/home/ole2001/PROGRAMS/ChIPseeqer-1.0/DATA/SOLEXA/100217_H3K4me3/CHIP";
my $format  = "eland";
my $files   = "*.matches";
my $submit  = 0;

GetOptions("files=s"   => \$files,
           "format=s"  => \$format,
	   "submit=s"  => \$submit,
	   "chipdir=s" => \$chipdir);


my $todo = undef;

my $a_ref_f = Sets::getFiles($files);

my $pwd  = `pwd`; $pwd =~ s/\n//;


foreach my $f (@$a_ref_f) {
  print "Processing $f\n";
  system("head -3 $f");

  my $pbs = PBS->new;
  $pbs->setPlatform("panda");
  $pbs->setWallTime("1:00:00");

  $pbs->addCmd("cd $pwd");
  $pbs->setScriptName("$f.script");    

  $todo = "perl ~/PERL_SCRIPTS/MyScanACE_matchesToIntervals.pl --motifmatches=$f --ext=200 > $f.extint200";
  $pbs->addCmd($todo);
  
  $todo = "ChIPseeqerGetReadDensityProfiles.bin -intervals $f.extint200 -chipdir $chipdir -format $format -fraglen 0 -uniquereads 1 > $f.extint200.H3K4";
  $pbs->addCmd($todo);
  
  $todo = "R --slave --args $f.extint200.H3K4 < /home/ole2001/PROGRAMS/ChIPseeqer/SCRIPTS/makeAvgReadProfile.R";
  $pbs->addCmd($todo);

  if ($submit == 1) {
    $pbs->submit();
  } else {
    $pbs->print();
  }

}
