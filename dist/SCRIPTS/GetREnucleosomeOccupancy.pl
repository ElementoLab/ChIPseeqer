#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args --lanes=FILE --motifs=FILE --chrdata=FILE --outfile=FILE \n";
}
my $lanes     = undef;
my $motifs    = undef;
my $chrdata   = undef;
my $datadir   = ".";
my $normalize = 0;
my $data      = undef;
my $verbose   = 0;
my $mapt      = 36;
my $mapdir    = undef;

#my $outfile   = undef;


GetOptions("lanes=s"     => \$lanes,
	   "normalize=s" => \$normalize,
           "motifs=s"    => \$motifs,
	   "data=s"      => \$data,
	   "mapdir=s"    => \$mapdir,
	   "mapt=s"      => \$mapt,
	#   "outfile=s"   => \$outfile,
	   "verbose=s"   => \$verbose,
	   "datadir=s"   => \$datadir,
	   "chrdata=s"   => \$chrdata);



my @lane_data = ();
my %M = ();
open IN, $lanes or die "Cannot open file $lanes\n";
while (my $l = <IN>) {
  chomp $l;
  my @ld = split /\t/, $l, -1;
  push @lane_data, \@ld;
  my $outfile = Sets::filename($motifs) . ".$ld[0]";
  my $prmfile = "$outfile.prm";
  my $todo = "$ENV{CHIPSEEQERDIR}/ChIPgetIntervalReadCounts -chipdir $datadir/$ld[0].fastq.sam_SPLIT -intervals $motifs -normalize $normalize -chrdata $chrdata -output max -format sam -outfile $outfile ";
  if (defined($mapdir)) {
    $todo .= " -mapdir $mapdir -mapt $mapt ";
  }
  if ($normalize == 1) {
    $todo .= " -normto 10000000 ";
  }
  system($todo) == 0 or die "Cannot exxecute $todo\n";

  if ($verbose == 1) {		    
    print "$todo\n";
  }	
  # get params
  open IN1, $prmfile or die "Cannot open file";
  my %T = ();
  while (my $l = <IN1>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    $T{$a[0]} = $a[1];
  }
  close IN1;

  # get 
  open IN1, $outfile or die "Cannot open file";
  while (my $l = <IN1>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    
    my $m = "$a[0]\t$a[1]\t$a[2]";
    $M{ $m }{ $ld[1] } = sprintf("%3.1f", $a[3]);
    
    if ($data eq "MAINE") {  # nucleosome
      if ($M{ $m }{ $ld[1] } <= $T{"q0.25"}) {
	$M{ $m }{ $ld[1] } .= "(F)";
      #} elsif ($M{ $m }{ $ld[1] } >= $T{"q0.75"}) {
	#$M{ $m }{ $ld[1] } .= "(N)";
      } else {
	$M{ $m }{ $ld[1] } .= "(-)";
      }
    }
      
    if ($data eq "FAIRE") {  # accessible
      if ($M{ $m }{ $ld[1] } > $T{"q0.25"}) {
	$M{ $m }{ $ld[1] } .= "(A)";
      } else { 
	$M{ $m }{ $ld[1] } .= "(-)";
      }
    }

  }
  close IN1;

}
close IN;

print "chr\tstart\tend";
foreach my $l (@lane_data) {
  print "\t$l->[1]"; 
}
print "\n";

foreach my $m (keys(%M)) {
  
  print "$m";
  foreach my $l (@lane_data) {
    print "\t$M{$m}{$l->[1]}"; 
  }
  print "\n";
}
