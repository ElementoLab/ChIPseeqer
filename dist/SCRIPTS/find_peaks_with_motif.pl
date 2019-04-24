#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}/SCRIPTS";
use lib "$ENV{CHIPSEEQERDIR}";

use Getopt::Long;
use Table;
use Sets;
use strict;

my $targets		= undef; # variable that stores the targets file
my $seed		= undef; # variable that stores the input seed
my $profiles	= undef; # variable that stores the profiles file
my $summary		= undef; # variable that stores the summary file
my $motif_expr	= undef; # variable that stores the motif expression
my $peak		= undef; # variable that stores the peak
my $outfile		= undef; # variable that stores the output file
my $fromfile	= undef; # variable that stores the original peaks file

my @seeds		= undef; # array that stores the seeds form the summary file
my @line		= undef; # array that stores a line from the summary file

# hashes
my %peaksHash	= ();

die "Please specify --targets=FILE --seed=STR \n" if (@ARGV == 0);


GetOptions("targets=s"	=> \$targets,
"seed=s"	=> \$seed,
"fromfile=s" => \$fromfile);

if(!defined($targets)) {
	die "Please specify --targets=FILE \n"
}
elsif(!defined($seed)) {
	die "Please specify --seed=STR \n"
}
else {
	$profiles	= "$targets.profiles";
	$summary	= "$targets.summary";
}

#
# Retrieve the expression motif from the summary file, given a seed
#

# open the seeds file 
open IN, "$summary";

# read seeds file line by line
while (my $l = <IN>) {
	chomp $l;
	@line = split /\t/, $l, -1;
	
	if($seed eq $line[8]) {
		$motif_expr = $line[0];
	}
}
close IN;

if(!defined($motif_expr)) {
	die "The seed you specified does not exist in these FIRE results. \n"
}

#
# Retrieve the peaks associated with the expression motif
#

# open profiles file 
open IN, "$profiles";
open OUTFILE, ">$targets.$seed.txt";

my $prevpeak = undef;
# read profiles file line by line
while (my $l = <IN>) {
	chomp $l;
	@line = split /\t/, $l, -1;
	
	# if the expression motif is found in the line
	if($motif_expr eq $line[0]) {
		$peak		= $line[1];
		next if $peak =~ /random$/;     # skip peak ending with random
		
		my @peakfield	= split /-/, $peak, -1;
		
		my $peakdist	= $peakfield[2] - $peakfield[1];		
		my $prct		= sprintf("%.1f", ($line[2]/$peakdist)*100);
		
		if(!defined($prevpeak)) {
			print OUTFILE "$peakfield[0]\t$peakfield[1]\t$peakfield[2]\t$line[2]:$prct";
		}
		else {
			if($peak eq $prevpeak) {
				print OUTFILE ",$line[2]:$prct";
			}
			else {
				print OUTFILE "\n";
				print OUTFILE "$peakfield[0]\t$peakfield[1]\t$peakfield[2]\t$line[2]:$prct";
			}
		}
	}
	$prevpeak = $peak;
}
close IN;
close OUTFILE;

print "Done ($targets.$seed.txt created).\n";

#
# Get peaks from a specific file
#
if (defined($fromfile)) {
	if(-e "$fromfile") {
		
		# open the original peaks file 
		open IN, "$fromfile";
		
		# store the input peaks in a hash
		while (my $l = <IN>) {
			chomp $l;
			my @ipeak	= split /\t/, $l, -1;
			my $tmp	= "$ipeak[0]\t$ipeak[1]\t$ipeak[2]";
			
			$peaksHash{$tmp} = $l;
		}
		close IN;
				
		# open the original peaks file 
		open IN, "$targets.$seed.txt";
		open OUTFILE2, ">$fromfile.$seed.txt";
		
		# store the input peaks in a hash
		while (my $l = <IN>) {
			chomp $l;
			my @mpeak	= split /\t/, $l, -1;
			my $tmp2	= "$mpeak[0]\t$mpeak[1]\t$mpeak[2]";
			
			if(exists($peaksHash{$tmp2})) {
				print OUTFILE2 "$peaksHash{$tmp2}\n";
			}			
		}
		close IN;
		close OUTFILE2;
		
		print "Done ($fromfile.$seed.txt created).\n";

	}
}
