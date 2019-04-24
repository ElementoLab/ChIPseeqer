#!/usr/bin/perl
use lib "$ENV{PERLMODULESDIR}";
use Sets;
use strict;

use Getopt::Long;

if (@ARGV == 0) {
	die "Args --promoters=FILE --peakfile=FILE [ --hicdata=FILE ] \n";
}
my $promoters			= "$ENV{CHIPSEEQERDIR}/DATA/hg18/refSeq.oneperTSS_u1000_d1000.txt";
my $peakfile			= undef;
my $verbose				= 0;
my $label				= "BINDING";
my $d0					= 5000;
my $peakheightcol		= undef;
my $hicdata				= undef;
my $minhicdist			= 2000;
my %HiC					= ();
my $haspeaksummit		= 1;

GetOptions("promoters=s"	=> \$promoters,
"verbose=s"					=> \$verbose,
"hicdata=s"					=> \$hicdata,
"label=s"					=> \$label,
"peakfile=s"				=> \$peakfile,
"d0=s"						=> \$d0,
"peakheightcol=s"			=> \$peakheightcol,
"haspeaksummit=s"			=> \$haspeaksummit);

if (defined($hicdata)) {
	
	open IN, $hicdata or die "Cannot open $hicdata\n";
	while (my $l = <IN>) {
		next if ($l =~ /^\#/);
			chomp $l;
		my @a = split /[\t\ +]/, $l, -1;
		$HiC{$a[0]}{"$a[2]\t$a[4]\t$a[5]"} = $a[7];    
	}
	close IN;
}

my %S = ();
open IN, $peakfile or die "Cannot open $peakfile\n";
while (my $l = <IN>) {
	chomp $l;
	my @a = split /\t/, $l, -1;
	$S{"$a[0]\t$a[1]\t$a[2]"} = $a[$peakheightcol-1];
}
close IN;


my $todo = "$ENV{CHIPSEEQERDIR}/CompareIntervals -peakfile1 $promoters -hasid1 1 -ext1 48000 -peakfile2 $peakfile -show_ov_int 1 -showsummit2 $haspeaksummit ";
my $txt  = `$todo`;

print "GENE\t$label\n";
my @a    = split /\n/, $txt;
foreach my $r (@a) {
		
	my @b = split /\t/, $r;
	my $d1 = ($b[2] + $b[3])/2;
	my $g = shift @b;
	
	shift @b;
	shift @b;
	shift @b;
	my $num = shift @b;
	
	my $sum = 0;
	
	for (my $i=0; $i<$num; $i++) {
		my $p = $b[$i];
		my @c = split /\-/, $p;
		my $d2 = undef;
		
		if($haspeaksummit) {
			$d2 = $c[3];
		}
		else {
			$d2 = ($c[2] + $c[1])/2;
		}
		
		my $d  = abs($d2 - $d1);
		my $pikid = "$c[0]\t$c[1]\t$c[2]";
		my $s  = $S{$pikid}; # score
		my $k = -1;
		my $score = undef;
		my $pen   = -1;
		if (!defined($hicdata)) {      
			$score = $s * exp(-$d/$d0);
		} else {
			if ($d <= 1000) { # low distance, use peak height
				$score  = $s;
			} elsif (defined($HiC{$g}{$pikid})) { # large distance
				$k = $HiC{$g}{$pikid};
				$pen = ( 1 / ( 1+(exp(-10*($k-1))) ) );
				#if ($k >= 2) {
				#  $pen = 1;
				#} else {
				#  $pen = 0;
				#}
				#if ($verbose == 1) {
				#  print "pen=$pen\n";
				#}
				$score = $s * $pen;
			} else {
				$pen = "undef";
				$score = 0;
			}
			
		}
		
		$sum += $score;
		if ($verbose == 1) {
			print " $s $d $pen $k $score\n";
		}
	}
	$sum = sprintf("%3.2f", $sum);
	print "$g\t$sum\n";
	
}
