#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args --matches=FILE --peakfile=FILE --chipdir=DIR --format=STR\n";
}

my $matches  = undef;
my $peakfile = undef;
my $todo     = undef;
my $format   = undef;
my $chipdir  = undef;
my $exec     = 1;
GetOptions("matches=s"   => \$matches,
	   "chipdir=s"   => \$chipdir,
	   "format=s"    => \$format,
	   "exec=s"      => \$exec,
           "peakfile=s"  => \$peakfile);
 
die "Please specify --peakfile\n" if !defined($peakfile);

#/home/ole2001/PROGRAMS/MYSCANACE/MyScanACE -jb /home/ole2001/PROGRAMS/ChIPseeqer-1.0/DATA/BULYK_MATRICES/Bcl6b_pwm_primary.txt -z /home/ole2001/PROGRAMS/SNPseeqer/REFDATA/hg18/wg.fa -c 1 -maxnummatches 100000 -intervals /home/ole2001/PROGRAMS/ChIPseeqer-1.0/DATA/SOLEXA/H3K4me3_NB_100526/H3K4me3NB_peaks_t15.txt 

$todo  = "perl $ENV{HOME}/PERL_SCRIPTS/MyScanACE_matchesToIntervals.pl --motifmatches=$matches --showall=1 > $matches.Int\n";
$todo .= "$ENV{CHIPSEEQERDIR}/ChIPgetIntervalReadCounts -intervals $matches.Int -chipdir $chipdir -format $format -normalize 0 -chrdata $ENV{CHIPSEEQERDIR}/DATA/hg18.chrdata > $matches.Int.H3K4\n";
$todo .= "$ENV{CHIPSEEQERDIR}/ChIPseeqerCons --peakfile=$matches.Int.H3K4 --consdir=$ENV{HOME}/PROGRAMS/ChIPseeqer-1.0/DATA/CONSERVATION/ --format=gzscores --category=placental --showalldata=1 --outfile=$matches.Int.H3K4.cons\n";
$todo .= "$ENV{CHIPSEEQERDIR}/CompareIntervals -peakfile1 $matches.Int.H3K4.cons       -peakfile2 $peakfile -showpeakdesc 1 > $matches.Int.H3K4.cons.peaks\n";
$todo .= "$ENV{CHIPSEEQERDIR}/CompareIntervals -peakfile1 $matches.Int.H3K4.cons.peaks -peakfile2 $ENV{CHIPSEEQERDIR}/DATA/CpGislands_annotation -showpeakdesc 1 >  $matches.Int.H3K4.cons.peaks.CGI\n";
$todo .= "perl $ENV{HOME}/PERL_SCRIPTS/FindMotifMatchesWithPattern.pl --matches=$matches.Int.H3K4.cons.peaks.CGI --pattern=\"..CG..\" --col=5 > $matches.Int.H3K4.cons.peaks.CGI.CG\n";

if ($exec == 1) {
  my $tmpfile = Sets::getTempFile("/tmp/gather");
  open OUT, ">$tmpfile";
  print OUT $todo;
  close OUT;
  system("bash $tmpfile");
} else {
  print "$todo";
}
