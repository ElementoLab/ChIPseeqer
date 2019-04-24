#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
#use lib "$ENV{HOME}/usr/lib/perl5/site_perl/5.8.8";
use Statistics::Test::WilcoxonRankSum;
use Sets;
use strict;
use Getopt::Long;

#use Statistics::TTest;

my $verbose   = 0;
my $exec      = 1;
my $stat      = 1;
my $format    = "eland";
my $chipdir   = "/home/ole2001/PROGRAMS/ChIPseeqer-1.0/DATA/SOLEXA/100217_H3K4me3/CHIP";
my $peakfile  = "/home/ole2001/PROGRAMS/ChIPseeqer-1.0/DATA/SOLEXA/100217_H3K4me3/h3k4me3_peaks.txt";
my $matches   = undef;
my $motif     = undef;
my $motiftype = "jb";
my $outdir    = ".";
my $todo      = undef;
my $genome    = "/home/ole2001/PROGRAMS/SNPseeqer/REFDATA/hg18/wg.fa";


if (@ARGV == 0) {
  die "ARgs: [ --matches=FILE | --motif=FILE --motiftype=STR (j or jb - default) --outdir=DIR (.) ] --peakfile=FILE --chipdir=FILE --format=STR --gendata=INT (0) \n";
}

GetOptions("matches=s"  => \$matches,
	   "gendata=s"  => \$exec,
	   "motif=s"    => \$motif,
	   "motiftype=s"=> \$motiftype,
	   "verbose=s"  => \$verbose,	
	   "format=s"   => \$format,
	   "outdir=s"   => \$outdir,
	   "genome=s"   => \$genome,
	   "chipdir=s"  => \$chipdir,
	   "peakfile=s" => \$peakfile); 	

if (! -e $genome) {
  die "$genome does not exist\n";
}

if (defined($motif) && defined($matches)) {
  die "Can't define both motif and matches\n";
}

if ($exec == 1) {

  if (defined($motif)) {    
    $matches = "$outdir/" . Sets::filename($motif) . ".matches";
    $todo = "/home/ole2001/PROGRAMS/MYSCANACE/MyScanACE -$motiftype $motif -z $genome -c 1 -maxnummatches 100000 -intervals $peakfile > $matches";
    print "$todo\n" if ($verbose == 1);
    system($todo) == 0 or die "Cannot exec $todo\n";
  }
  
  $todo    = "perl /home/ole2001/PERL_SCRIPTS/MyScanACE_matchesToIntervals.pl --motifmatches=$matches --ext=0 --showaff=1  > $matches.Int";
  print "$todo\n" if ($verbose == 1);
  system($todo) == 0 or die "Cannot exec $todo\n";
  
  $todo       = "/home/ole2001/PROGRAMS/ChIPseeqer/ChIPgetIntervalReadCounts -intervals $matches.Int -chipdir $chipdir -format $format -normalize 0 -chrdata /home/ole2001/PROGRAMS/ChIPseeqer/DATA/hg18.chrdata > $matches.Int.H3K4";
  print "$todo\n" if ($verbose == 1);
  system($todo) == 0 or die "Cannot exec $todo\n";
  
  $todo       = "perl /home/ole2001/PROGRAMS/ChIPseeqer/SCRIPTS/GenerateRandomMotifOccurrencesInPeaks.pl --peakfile=$peakfile --motifmatches=$matches.Int > $matches.Int.random";
  print "$todo\n" if ($verbose == 1);
  system($todo) == 0 or die "Cannot exec $todo\n";
  
  $todo       = "/home/ole2001/PROGRAMS/ChIPseeqer/ChIPgetIntervalReadCounts -intervals $matches.Int.random -chipdir $chipdir -format $format -normalize 0 -chrdata /home/ole2001/PROGRAMS/ChIPseeqer/DATA/hg18.chrdata > $matches.Int.random.H3K4";
  print "$todo\n" if ($verbose == 1);
  system($todo) == 0 or die "Cannot exec $todo\n";
}

my @th = ();
open IN, "$matches.Int.H3K4" or die "Cannot open H3K4 file\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @th, $a[4];
}
close IN;

my @rh = ();
open IN, "$matches.Int.random.H3K4" or die "Cannot open rand file\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @rh, $a[3];
}
close IN;

my $m1 = Sets::median(\@rh);
my $m2 = Sets::median(\@th);
my $ra = $m1/$m2;
print "$matches\t$m1\t$m2\t$ra";
if ($stat== 1) {
  my $wilcox_test = Statistics::Test::WilcoxonRankSum->new;
  $wilcox_test->load_data(\@rh, \@th);  
  my $prob = $wilcox_test->probability();
  my $pf = sprintf '%5.3e', $prob; 
  print "\t$pf";
  #my $ttest = new Statistics::TTest;
  #$ttest->load_data(\@rh,\@th);
  #$ttest->output_t_test();
}
print "\n";

#ole2001@panda BULYK_MATRICES $ columns.pl 3 <  Bcl6b_pwm_primary.txt.matches.Int.random.H3K4 | getMedian.pl -
#15.273
#ole2001@panda BULYK_MATRICES $ columns.pl 4 <  Bcl6b_pwm_primary.txt.matches.Int.H3K4 | getMedian.pl -
#12.636
