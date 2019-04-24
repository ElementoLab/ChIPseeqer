#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Getopt::Long;

# variables to store the arguments values
my $matrixfile		= undef;
my $sortby			= undef;
my $suf				= "sorted";
my $verbose			= 0;

# handling lack of arguments
if (@ARGV == 0) {
	die "Usage: sort_matrix.pl --matrixfile=FILE --sortby=FILE --suf=STR --verbose=INT \n";
}

# processing command line options
GetOptions("matrixfile=s" => \$matrixfile,
"sortby=s"		=> \$sortby,
"suf=s"			=> \$suf,
"verbose=i"		=> \$verbose );

if (!defined($suf)) {
	die("Must provide --suf=STR\n");
}


#
# Replace first tab with space
#
my $todo="perl -pi -e 's/\\t/ /' $matrixfile";
if ($verbose == 1) {
	print "$todo\n";
}
system($todo) == 0 or die "Cannot execute $todo.\n"; 


#
# Grep line from matrix file
#
$todo="for l in `cat $sortby`; do grep \"^\$l \" $matrixfile ; done > $matrixfile.$suf";
if ($verbose == 1) {
	print "$todo\n";
}
print "Sorting $matrixfile by $sortby...";
system($todo) == 0 or die "Cannot execute $todo.\n"; 
print "Done.\n";


#
# Replace first space with tab
#
if (-e "$matrixfile.$suf") {
	$todo= "perl -pi -e 's/ /\\t/' $matrixfile";
	if ($verbose == 1) {
		print "$todo\n";
	}
	system($todo) == 0 or die "Cannot execute $todo.\n";
}