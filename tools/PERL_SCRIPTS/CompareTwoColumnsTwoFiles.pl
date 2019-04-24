#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
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

my $file1     = undef;
my $col1      = undef;
my $file2     = undef;
my $col2      = undef;

if (@ARGV == 0) {
  die "ARgs:  --file1=FILE --col1=INT --file2=FILE --col2=INT \n";
}

GetOptions("file1=s"    => \$file1,
	   "file2=s"    => \$file2,
	   "col1=s"     => \$col1,
	   "col2=s"     => \$col2);

my @th = ();
open IN, "$file1" or die "Cannot open $file1\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @th, $a[$col1];
}
close IN;

my @rh = ();
open IN, "$file2" or die "Cannot open $file2\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @rh, $a[$col2];
}
close IN;

my $m1 = Sets::median(\@rh);
my $m2 = Sets::median(\@th);
my $ra = $m1/$m2;
print "$file1\t$file2\t$m1\t$m2\t$ra";

my $wilcox_test = Statistics::Test::WilcoxonRankSum->new;
$wilcox_test->load_data(\@rh, \@th);  
my $prob = $wilcox_test->probability();
my $pf = sprintf '%5.3e', $prob; 
print "\t$pf";

print "\n";

