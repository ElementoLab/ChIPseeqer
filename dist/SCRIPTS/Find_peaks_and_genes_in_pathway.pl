#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";
use lib "$ENV{PAGEDIR}";

use Getopt::Long;
use Table;
use Sets;
use strict;

my $expfile		= undef;
my $bin			= undef;
my $pathway		= undef;
my $species		= undef;
my $onecatfile	= undef;
my $exptype		= undef;
my $verbose		= undef;
my $pref		= undef;
my $overlapfile = undef;

# hashes
my %NMHash		= ();

# handling of missing arguments
if (@ARGV == 0) {
	die "Usage: perl find_peaks_and_genes_in_pathway.pl --expfile=FILE --bin=INT --pathway=STR --species=STR (e.g. human_go_orf) --prefix=STR\n";
}

GetOptions("expfile=s"  => \$expfile,
"pathway=s"  => \$pathway,
"bin=s"      => \$bin,
"exptype=s"  => \$exptype,
"onecatfile=s" => \$onecatfile,
"species=s"  => \$species,
"verbose=s" => \$verbose,
"prefix=s"	=> \$pref,
"overlapfile=s" => \$overlapfile);

if (!defined($pref)) {
	die("Must provide --prefix=STR\n");
}

if (!defined($pref)) {
	die("Must provide --overlapfile=STR\n");
}

my $todo = "perl $ENV{PAGEDIR}/SCRIPTS/find_genes_in_bin_and_pathway.pl --expfile=$expfile --bin=$bin --pathway=$pathway --species=$species > $pref.genes.$pathway.txt";
if ($verbose == 1){ 
	print "$todo\n";
}
print "Looking for genes in bin $bin and pathway $pathway ...";
system($todo) == 0 or die "Cannot exec $ENV{PAGEDIR}/SCRIPTS/find_genes_in_bin_and_pathway.pl\n";

if (-e "$pref.genes.$pathway.txt") {
	print "Done ($pref.genes.$pathway.txt created).\n";	
}

print "Looking for peaks in bin $bin and pathway $pathway ...";

# fill in the NMs hash from the overlapfile
open IN, "$overlapfile";

my $cnt = 0;
while (my $l = <IN>) {
	chomp $l;
	my @a = split /\t/, $l, -1;
	
	my $peak	= join("\t", ($a[0], $a[1], $a[2]));
	
	
	for(my $i=3; $i<$#a+1; $i++) {
		my @b = split /\t/, $a[$i];
		
		$NMHash{$b[0]} = $peak;
	}	
}
close IN;

#for my $key ( keys %NMHash ) {
#	my $value = $NMHash{$key};
#	print "$key => $value\n";
#}

# looking for the peaks for each NM
open IN, "$pref.genes.$pathway.txt";
open PEAKSFILE, ">$pref.peaks.$pathway.txt";

my $cnt = 0;
while (my $l = <IN>) {
	chomp $l;
	
	#split the line
	my @a = split /\t/, $l, -1;
	
	my $NMs = $a[1];
	
	#split the NMs field
	my @b = split /\//, $NMs, -1;
	
	foreach my $j(@b) {
		if(exists($NMHash{$j})) {
			print PEAKSFILE "$NMHash{$j}\t$j\n";
		}
	}
	
}
close IN;
close PEAKSFILE;

print "Done ($pref.peaks.$pathway.txt created).\n";	

