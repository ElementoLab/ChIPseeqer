#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";
use Sets;
use strict;
use Getopt::Long;

my $file1	= undef;
my $file2	= undef;
my $type	= "union"; # Can be union or unique 

GetOptions
(	"peakfile1=s"   => \$file1,
"peakfile2=s"   => \$file2,
"type=s"		=> \$type);

if (! -e $file1) {
	die "$file1 cannot be opened\n";
} 
if (! -e $file2) {
	die "$file2 cannot be opened\n";
}


srand();

my $ran      = int(rand(5782763));  
my $tmpfile1 = "tmpfile1.$ran";

if($type eq "union") {
	
	my $todo = "$ENV{CHIPSEEQERDIR}/CompareIntervals -peakfile1 $file1 -peakfile2 $file2 -showunion 1 -show_ov_int 1 > $tmpfile1";
	system($todo) == 0 or die "Cannot exec $todo\n";
	
	open IN, $tmpfile1;
	
	while (my $l = <IN>) {
		chomp $l; 
		my @c = (split /\t/, $l); 
		
		if ($c[3] == 0) {
			print "$c[0]\t$c[1]\t$c[2]\n";
		} elsif ($c[3] >= 1) {
			my $last	= (split /\t/, $l)[-1]; 
			my @d		= (split /-/, $last);
			print "$d[0]\t$d[1]\t$d[2]\n";
			
		}
	}
	
	close IN;
}
elsif($type eq "unique") {
	my $todo = "$ENV{CHIPSEEQERDIR}/CompareIntervals -peakfile1 $file1 -peakfile2 $file2 > $tmpfile1";
	system($todo) == 0 or die "Cannot exec $todo\n";
	
	open IN, $tmpfile1;
	
	while (my $l = <IN>) {
		chomp $l; 
		my @c = (split /\t/, $l); 
		
		if ($c[3] == 0) {
			print "$c[0]\t$c[1]\t$c[2]\n";
		}
	}
	
	close IN;
}

#the inverse compare intervals
srand (time);     

my $ran      = int(rand(5782763)); # generate random suffix   
my $tmpfile2 = "tmpfile2.$ran";   

my $todo = "$ENV{CHIPSEEQERDIR}/CompareIntervals -peakfile1 $file2 -peakfile2 $file1 > $tmpfile2";
system($todo) == 0 or die "Cannot exec $todo\n";

open IN, $tmpfile2;

while (my $l = <IN>) {
	chomp $l; 
	my @c = (split /\t/, $l); 
	if ($c[3] == 0) {
		print "$c[0]\t$c[1]\t$c[2]\n";
	}	
}

close IN;

unlink $tmpfile1;
unlink $tmpfile2;

