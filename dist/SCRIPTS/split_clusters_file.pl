#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use FileHandle;
use Getopt::Long;
use strict;

# variables to store the arguments values
my $clustersfile		= undef;
my $clustersnumber		= undef;

my %FH		= ();
my $i		= 0;
my $todo	= undef;

# handling lack of arguments
if (@ARGV == 0) {
	die "Usage: split_clusters_file.pl --clustersfile=FILE --clustersnumber=INT \n";
}

# processing command line options
GetOptions("clustersfile=s" => \$clustersfile,
"clustersnumber=i"			=> \$clustersnumber);

if (!defined($clustersnumber)) {
	die("Must provide --clustersnumber=INT\n");
}

# open clustersfile
open IN, $clustersfile;

# for each line
while (my $l = <IN>) {
	chomp $l;
	
	#split line
	my @a = split /[\t]/, $l, -1;
	
	my $name			= $a[0];
	my $clusterindex	= $a[1];
	
	#create a new clusters file
	my $fh = undef;
	
	if (!defined($FH{$clusterindex})) {
		$FH{$clusterindex} = new IO::File "> cluster.$clusterindex.txt";
	}
	
	$fh = $FH{$clusterindex};
	
	#write in reads file
	print $fh "$name\t$clusterindex\n";
}
close IN;

# put all cluster files in folder
$todo = "mkdir CLUSTERS";
system($todo) == 0 or die "Cannot execute $todo.\n";

$todo = "for l in `ls cluster*`; do mv \"\$l\" CLUSTERS ; done";
system($todo) == 0 or die "Cannot execute $todo.\n";