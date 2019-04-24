#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Fasta;
use Getopt::Long;
use strict;

# variables to store options
my $matrix1	= undef;
my $matrix2	= undef;

# other variables
my $verbose  = 0;
my @TF1		 = ();
my @TF2		 = ();
my %LABELS1  = ();
my %LABELS2  = ();

# handling of missing arguments
if (@ARGV == 0) {
	die "Usage: CompareJaccardIndexMatrices.pl --matrix1=FILE --matrix2=FILE\n";
}

# handling given options
GetOptions("matrix1=s"		=> \$matrix1,
"matrix2=s"		=> \$matrix2,
"verbose=s"		=> \$verbose);

if (!defined($matrix1)) {
	die("Must provide --matrix1=FILE\n");
}

if (!defined($matrix2)) {
	die("Must provide --matrix2=FILE\n");
}

#
# Make new files
#

# Open files
open(INM1, "$matrix1") or die "Can't open file $matrix1.";
open(INM2, "$matrix2") or die "Can't open file $matrix2.";

# skip 1st lines 
my $l1 = <INM1>; chomp $l1;
my $l2 = <INM2>; chomp $l2;

# store the matrix1 values in hash
while (my $line1 = <INM1>) {
	
	chomp $line1;
	
	my @b = split /\t/, $line1, -1;
	
	my $label = shift @b;

	my $values = join("\t", @b); 
	
	if($verbose==1) {
		print "Opening file $matrix1.\n";
	}
	
	push @TF1, $values;
	$LABELS1{$values} = $label;

}

# store the matrix2 values in hash
while (my $line2 = <INM2>) {
	
	chomp $line2;
	
	my @b = split /\t/, $line2, -1;
	
	my $label = shift @b;
	
	my $values = join("\t", @b); 
	
	if($verbose==1) {
		print "Opening file $matrix2.\n";
	}
	
	push @TF2, $values;
	$LABELS2{$values} = $label;
}

# close files
close INM1;
close INM2;

# Print the labels for the TFs
foreach my $f (@TF1) {
	print "\t$LABELS1{$f}";
}
print "\n";

my $idx = 0;

# for each line of matrix1
foreach my $tf (@TF1) { 
	
	my @tf1 = split /\t/, $tf, -1;
	my @tf2 = split /\t/, $TF2[$idx], -1;

	my $count = 0;
	
	print "$LABELS1{$tf}";
	
	foreach my $f(@tf1) {

		my $ratio = ($tf1[$count]+1) / ($tf2[$count]+1);
		my $logratio = log($ratio);

		#print "$LABELS1{$tf}\t$tf1[$count]\t$tf2[$count]\t$ratio\t$logratio\n";
		print "\t$logratio";

		$count++;
	}
	print "\n";
	
	$idx++;
}

print "\n";

