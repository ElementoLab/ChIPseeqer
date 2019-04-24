#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Getopt::Long;

# variables to store the arguments values
my $matrixfile		= undef;
my $sortby			= undef;
my $suf				= "sorted";
my $targets			= undef;
my $lenu			= 2000;
my $lend			= 2000;
my $generegion		= "TSS";	
my $prefix			= undef;
my $chipdir			= undef;
my $verbose			= 0;
my $matrixsuf		= "density";
my $cdt				= "cdt";

# handling lack of arguments
if (@ARGV == 0) {
	die "Usage: make_sort_cdt_matrix.pl --targets=FILE --lenu=INT --lend=INT --generegion=STR --prefix=STR --chipdir=STR --sortby=FILE --suf=STR --verbose=INT \n";
}

# processing command line options
GetOptions("sortby=s"		=> \$sortby,
"suf=s"			=> \$suf,
"targets=s"		=> \$targets,
"lenu=i"		=> \$lenu,
"lend=s"		=> \$lend,
"generegion=s"	=> \$generegion,
"prefix=s"		=> \$prefix,
"chipdir=s"		=> \$chipdir,
"format=s"		=> \$format,
"verbose=i"		=> \$verbose );

if (!defined($suf)) {
	die("Must provide --suf=STR\n");
}


#
# Create DensityMatrix file
#
my $ todo = "perl ChIPseeqer2DensityMatrix --targets=$targets --lenu=$lenu --lend=$lend --generegion=$generegion --prefix=$prefix --chipdir=$chipdir";
if ($verbose == 1) {
	print "$todo\n";
}
system($todo) == 0 or die "Cannot execute $todo.\n"; 


#
# Sort DensityMatrix file by another file
#
$todo = "perl SCRIPTS/sort_matrix.pl --matrixfile=$prefix.$matrixsuf --sortby=$sortby";
if ($verbose == 1) {
	print "$todo\n";
}
system($todo) == 0 or die "Cannot execute $todo.\n"; 


#
# Replace first space with tab
#
if (-e "$prefix.$matrixsuf.$suf") {
	$todo= "perl SCRIPTS/expression_matrix_to_CDT.pl $prefix.$matrixsuf.$suf > $prefix.$matrixsuf.$suf.$cdt";
	if ($verbose == 1) {
		print "$todo\n";
	}
	print "Creating cdt file... ";
	system($todo) == 0 or die "Cannot execute $todo.\n";
	if (-e "$prefix.$matrixsuf.$suf.$cdt") {
		print "Done ($prefix.$matrixsuf.$suf.$cdt created).\n";	
	}
}