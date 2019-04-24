#!/usr/bin/perl
use lib "$ENV{PERLMODULESDIR}";
use Sets;
use strict;

use Getopt::Long;

if (@ARGV == 0) {
	die "Args --genefile=FILE --geneparts=STR [ --showORF=INT --showprofile=INT ]\n";
}

my $genefile	= undef;
my $geneparts	= "P"; # can be "P,I", or any combination or P, E, I, I1, I2, DW, D, IG
my $showORF		= 0;
my $showprofile = 0;
my $outfile		= undef;

GetOptions(
"genefile=s"    => \$genefile,
"geneparts=s"	=> \$geneparts,
"showprofile=s" => \$showprofile,
"showORF=s"     => \$showORF,
"outfile=s"		=> \$outfile);

if (!defined($outfile)) {
	die("Must provide --outfile=STR\n");
}

my @a_p = split /\,/, $geneparts;

open IN, $genefile or die "Cannot open $genefile\n";
open OUTFILE, ">$outfile";

if ($showprofile == 1) {
	print OUTFILE "GENE\tTGT\n";
}

my $l = <IN>; chomp $l;
my @h = split /\t/, $l;
while (my $l = <IN>) {
	chomp $l;
	my @a = split /\t/, $l, -1;
	
	my $num = 0;
	for (my $i=1; $i<@a; $i++) {
		if (Sets::in_array($h[$i], @a_p)) {
			if ($a[$i] > 0) {
				$num ++;
			}
		}
	}
	
	if ($num > 0) {
		if($showORF) {
			print OUTFILE "$a[9]";
		}
		else{
			print OUTFILE "$a[0]";
		}
		if ($showprofile == 1) {
			print OUTFILE "\t1";
		}
		print OUTFILE "\n";
	} else {
		if ($showprofile == 1) {
			if($showORF) {
				print OUTFILE "$a[9]\t0\n";
			} else {
				print OUTFILE "$a[0]\t0\n";
			}
		}		
	}
	
}
close IN;

