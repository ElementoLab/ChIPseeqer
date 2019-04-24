#!/usr/bin/perl

#
# This script is used to create the index.txt and names.txt files needed 
# for iPAGE, for the Reactome pathways (ReactomePathways.gmt)
#

use Getopt::Long;
use strict;

my $database		= undef;
my $verbose			= undef;

my %pathways_idx	= ();
my %genes_idx		= ();

# handling of missing arguments
if (@ARGV == 0) {
	die "Usage: reactome_get_index_and_names.pl --database=FILE \n";
}

GetOptions("database=s"  => \$database);

if (!defined($database)) {
	dir("Please define --database=FILE.\n");
}

open INDEX_FILE, ">$database\_index.txt";
open NAMES_FILE, ">$database\_names.txt";


open IN, "$database";

my $i = 0; # counter

while (my $l = <IN>) {
	chomp $l;
	my @a = split /\t/, $l, -1;

	my $pathway = $a[0];
	
	print NAMES_FILE "REACT_$i\t$pathway\n";
	$pathways_idx{$pathway}= "REACT_$i";
	
	for(my $j=2; $j<$#a+1; $j++) {
		my @b = split /\t/, $a[$j];
		
		my $genename	= $b[0];
		
		if(exists($genes_idx{$genename})) {
			$genes_idx{$genename} = join("\t", ($genes_idx{$genename}, $pathways_idx{$pathway}));
		}
		else {
			$genes_idx{$genename} = $pathways_idx{$pathway};
		}
	}
		
	$i++;
}
close IN;

#for my $key ( keys %pathways_idx ) {
#	print  "$key\t$pathways_idx{$key}\n";
#}

for my $key ( keys %genes_idx ) {
	print INDEX_FILE "$key\t$genes_idx{$key}\n";
}

close INDEX_FILE;
close NAMES_FILE;
