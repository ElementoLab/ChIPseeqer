#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Getopt::Long;

# variables to store the arguments values
my $matrixfile		= undef;
my $middleblack		= 1;
my $normalizerow	= 1;
my $bgcolor			= undef;	
my $mincolor		= undef;
my $maxcolor		= undef;
my $missingcolor	= undef;
my $size			= "4:1";	
my $contrast		= 10;
my $numcolors		= 64;
my $suf				= undef;
my $verbose			= 1;

# handling lack of arguments
if (@ARGV == 0) {
	die "Usage: matrix2pngheatmap.pl --matrixfile=FILE --middleblack=INT --normalizerow=INT --bgcolor=STR --mincolor=STR --maxcolor=STR --missingcolor=STR --size=STR --contrast=INT --numcolors=INT --suf=STR --verbose=INT \n";
}

# processing command line options
GetOptions("matrixfile=s" => \$matrixfile,
"middleblack=i"		=> \$middleblack,
"normalizerow=i"	=> \$normalizerow,
"bgcolor=s"			=> \$bgcolor,
"mincolor=s"		=> \$mincolor,
"maxcolor=s"		=> \$maxcolor,
"missingcolor=s"	=> \$missingcolor,
"size=s"			=> \$size,
"contrast=i"		=> \$contrast,
"numcolors=i"		=> \$numcolors,
"suf=s"				=> \$suf,
"verbose=i"			=> \$verbose );

if (!defined($suf)) {
	die("Must provide --suf=STR\n");
}

#
# Run matrix2png
#
my $todo = "matrix2png -data $matrixfile -b -z -bkgcolor $bgcolor -mincolor $mincolor -maxcolor $maxcolor -missingcolor $missingcolor -size $size -con $contrast"  ;
if ($middleblack == 1) {
	$todo .= " -b ";
}

if ($normalizerow == 1) {
	$todo .= " -z ";
}

$todo .= " > $suf";

if ($verbose == 1) {
	print "$todo\n";
}

print "Creating heatmap $suf...";
system($todo) == 0 or die "Cannot execute matrix2png program.\n"; 
print "Done.\n";