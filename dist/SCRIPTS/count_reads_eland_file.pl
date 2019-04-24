#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Sets;
use FileHandle;
use strict;

#handling lack of arguments
if (@ARGV == 0) {
  die "Args: readfiles\n";
}

#number of unique reads
my $unum = 0;

#number of multiple (repeated) reads
my $rnum = 0;

#number of total reads
my $tnum = 0;

#ratio of multiple reads to total reads
my $ratio =0;


#for every file
foreach my $f (@ARGV) {
  
  print STDERR "Opening $f\n";

  #open file
  open IN, $f;
  
  #for every line
  while (my $l = <IN>) {
    chomp $l;
    
	#split the line  
	my @a = split /[\t]/, $l, -1;
	
	#ignore the not matched reads
    if (($a[2] eq 'NM') || ($a[2] eq 'QC') || ($a[2] eq 'RM')) {
		next;
	}
	#count a unique read
	elsif (($a[2] eq 'U0') || ($a[2] eq 'U1') || ($a[2] eq 'U2')) {
		$unum++;
	}
	#count a multiple read
	elsif (($a[2] eq 'R0') || ($a[2] eq 'R1') || ($a[2] eq 'R2')) {
		$rnum++;
	}
  }
	
	#calculate the number of total reads fount (unique and multiple)
	$tnum = $unum + $rnum;

	#calculate the ratio
	if($tnum != 0) {
		$ratio = $rnum/$tnum;
	}
	
	print "Found $unum unique reads	\n";
	print "Found $rnum multiple reads \n";
	print "Found $tnum total reads \n";
	print "Non-unique (multiple) reads ratio is $ratio";
	
	print "\n";
	close IN;
	
}