#!/usr/bin/perl
use lib "$ENV{CHIPSEEQERDIR}";

use Sets;
use Table;
use Getopt::Long;
use strict;

my $peakfile = undef;
my $w        = 2000;
my $force    = 0;
my $nosummit = 0;

if (@ARGV == 0) {
	die "Args: --peakfile=FILE --w=INT\n";
}

GetOptions("peakfile=s" => \$peakfile,
"force=s"    => \$force,
"w=s"        => \$w,
"nosummit=s" => \$nosummit);

my $outfile = "$peakfile.centered$w";
if (($force == 0) && (-e $outfile)) {
	print "Warning: $outfile exists. Overwrite ?\n";
	<STDIN>;
}

open OUT, ">$outfile";
my $ta = Table->new;
$ta->loadFile($peakfile);

if($nosummit == 0) {
	$ta->sortbycol(4);
	my $a_ref = $ta->getArray();
	
	foreach my $r (@$a_ref) {
		
		my $i = $r->[5] - $w/2;
		my $j = $r->[5] + $w/2;
		next if ($i < 0);
		
		print OUT "$r->[0]\t$i\t$j\n";
	}
	
}
else {
	my $a_ref = $ta->getArray();

	foreach my $r (@$a_ref) {
		
		# midpoint = $r->[2] - $r->[1] /2
		# i = midpoint - w/2
		# j = midpoint w/2 
		
		#print OUT "$r->[0]\t$i\t$j\n";
		#print "$r->[0]\t$i\t$j\n";
	}

}
close OUT;

print "$outfile written\n";

