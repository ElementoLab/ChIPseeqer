#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use FileHandle;
use strict;

#handling lack of arguments
if (@ARGV == 0) {
	die "Args: file (in bed format)\n";
}

#open bed file
open IN, $ARGV[0];

#create new file
my $fh = new IO::File ">$ARGV[0].eland";

#for every line in the bed file
while (my $l = <IN>) {
	chomp $l;
	
	#split the line
	my @a		= split /[\ \t]/, $l, -1;
	
	#get the sequence name
	my $seqname = $a[3];
	
	#set the sequence
	my $seq		= "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
	
	#get the chromosome
	my $chr		= $a[0];
	
	#get the chromosome start position
	my $pos		= $a[1];
	
	#get the strand
	my $str		= undef;
		
	if ($a[5] eq "+") {
		$str = "F";
	}
	elsif ($a[5] eq "-") {
		$str = "R";
	}
	
	#write in eland file
	print $fh "$seqname\t$seq\tU0\t1\t0\t0\t$chr\t$pos\t$str\n"; 
	
}
close IN;