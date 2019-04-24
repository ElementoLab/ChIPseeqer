#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Fasta;
use Getopt::Long;
use strict;

# variables to store input files
my $input1		= undef;
my $input2		= undef;

# variables to store path to output files
my $output		= undef;

# variables to store lines
my $input1Line	= undef;
my $input2Line	= undef;

# hashes
my %input1		= ();
my %input2		= ();

# variables to store the numbers of lines
my $uniqcnt1	= 0;
my $uniqcnt2	= 0;
my $commcnt		= 0;
my $keepUniqRec	= 0; #if set to 1 eliminates double records

# variable to store other options
my $guiVersion	= 0;

# handling of missing arguments
if (@ARGV == 0) {
	die "Usage: CompareGenes.pl --input1=FILE --input2=FILE --output=PATH\n";
}

# handling given options
GetOptions("input1=s"	=> \$input1,
"input2=s"				=> \$input2,
"output=s"				=> \$output,
"keepUniqRec=s"			=> \$keepUniqRec,
"guiVersion=s"			=> \$guiVersion);

if (!defined($input1)) {
	die("Must provide --input1=FILE\n");
}

if (!defined($input2)) {
	die("Must provide --input2=FILE\n");
}

if (!defined($output)) {
	die("Must provide --output=PATH\n");
}

#
# open input1/input2 files and store in hashes
#
open(IN, "$input1") or die "Can't open file $input1.";

while ($input1Line = <IN>) {
	
	chomp $input1Line;
	
	# fill in input hash
	if($keepUniqRec == 1) {
		if(!exists($input1{$input1Line})) {
			$input1{$input1Line} += 1;
		}
	} else {
		$input1{$input1Line} += 1;
	}
}
close IN;

open(IN, "$input2") or die "Can't open file $input2.";

while ($input2Line = <IN>) {
	
	chomp $input2Line;
	
	# fill in input hash
	if($keepUniqRec == 1) {
		if(!exists($input2{$input2Line})) {
			$input2{$input2Line} += 1;
		}
	} else {
		$input2{$input2Line} += 1;
	}
}
close IN;

# make new files
open COMMONS, ">$output\Commons.txt";
open ALLUNIQUE, ">$output\AllUnique.txt";
open UNIQUE1, ">$output\uunique1.txt";
open UNIQUE2, ">$output\uunique2.txt";

#
# look for genes from input1
#
if($keepUniqRec == 1) {
	
	for my $key ( keys %input1 ) {
		
		if(exists($input2{$key})) {
			print COMMONS "$key\n";
		}
		else {
			print UNIQUE1 "$key\n";
			print ALLUNIQUE "$key\n";
		}
	}
	
}
else {
	open(IN, "$input1") or die "Can't open file $input1.";
	
	while ($input1Line = <IN>) {
		
		chomp $input1Line;
		if(exists($input2{$input1Line})) {
			print COMMONS "$input1Line\n";
		}
		else {
			print UNIQUE1 "$input1Line\n";
			print ALLUNIQUE "$input1Line\n";
		}
	}
	close IN;
}

close COMMONS;
close UNIQUE1;

#
# look for genes from input2
#
if($keepUniqRec == 1) {
	
	for my $key ( keys %input2 ) {
		
		if(!exists($input1{$key})) {
			print UNIQUE2 "$key\n";
			print ALLUNIQUE "$key\n";
		}
	}	
}
else {
	open(IN, "$input2") or die "Can't open file $input2.";
	
	while ($input2Line = <IN>) {
		
		chomp $input2Line;
		if(!exists($input1{$input2Line})) {
			print UNIQUE2 "$input2Line\n";
			print ALLUNIQUE "$input2Line\n";
		}
	}
	close IN;
}

close ALLUNIQUE;
close UNIQUE2;


#count number of unique and common genes
open IN, "$output\uunique1.txt";
while (my $l = <IN>) {
	chomp $l;
	my $line = split /\n/;
	$uniqcnt1++;
}
close IN;

open IN, "$output\uunique2.txt";
while (my $l = <IN>) {
	chomp $l;
	my $line = split /\n/;
	$uniqcnt2++;
}
close IN;

open IN, "$output\Commons.txt";
while (my $l = <IN>) {
	chomp $l;
	my $line = split /\n/;
	$commcnt++;
}
close IN;

if($guiVersion == 0) {
	print "Found:\t $uniqcnt1 unique genes in file $input1\n";
	print "\t $uniqcnt2 unique genes in file $input2\n";
	print "\t $commcnt common genes in files $input1 and $input2\n";
}