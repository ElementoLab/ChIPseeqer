#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";
use Sets;
use strict;
use Getopt::Long;

# variables to store options
my $file1 = undef;
my $file2 = undef;

# handling given options
GetOptions
(	"peakfile1=s"   => \$file1,
	"peakfile2=s"   => \$file2);

if (! -e $file1) {
	die "$file1 cannot be opened\n";
} 
if (! -e $file2) {
	die "$file2 cannot be opened\n";
}

#
# Getting the intersection of the two peak files
#
srand();

# Run CompareIntervals one way
my $ran      = int(rand(5782763));  
my $tmpfile3 = "tmpfile3.$ran";

my $I1 = "$ENV{CHIPSEEQERDIR}/CompareIntervals -peakfile1 $file1 -peakfile2 $file2 -output peaklist -ovtype AND > $tmpfile3";
system($I1) == 0 or die "Cannot exec $I1\n";
	
open IN, $tmpfile3;

my $numlines1 = 0;

while (my $l = <IN>){
	chomp $l; 
	 $numlines1++;
#	 print "$numlines1\n"; 
}
close IN; #first intersection coefficient

# Run CompareIntervals the other way too
my $ran      = int(rand(5782763));  
my $tmpfile4 = "tmpfile4.$ran";

my $I2 = "$ENV{CHIPSEEQERDIR}/CompareIntervals -peakfile1 $file2 -peakfile2 $file1 -output peaklist -ovtype AND > $tmpfile4";
system($I2) == 0 or die "Cannot exec $I2\n";	
	
open IN, $tmpfile4;
	my $numlines2=0; 

while (my $l = <IN>){
	chomp $l;
	$numlines2++;
#	print "$numlines2\n";
}
close IN; #second intersection coefficient


#
# Getting the union of the two peak files
#
my $ran      = int(rand(5782763));  
my $tmpfile5 = "tmpfile5.$ran";

my $union = "perl $ENV{CHIPSEEQERDIR}/CompareIntervalsOutputMergedList.pl -peakfile1 $file1 -peakfile2 $file2 > $tmpfile5";
system($union) == 0 or die "Cannot exec $union\n";	

open IN, $tmpfile5;
my $numlines3=0;

while (my $l = <IN>){
	chomp $l;
	$numlines3++;
#	print "$numlines3\n"
}
close IN;	

#
# Compute the Jaccard Index
#
my $J1= $numlines1/$numlines3;
my $J2= $numlines2/$numlines3;
#print "$numlines3\n";
#print "$numlines1\n";
#print "$numlines2\n";

#the average of the two J coefficients

my $Javg= ($J1+$J2)/2;
print sprintf("%3.2f\t", $Javg);

unlink $tmpfile3;
unlink $tmpfile4;
unlink $tmpfile5;