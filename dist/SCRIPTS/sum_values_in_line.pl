#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$ENV{CHIPSEEQERDIR}";

open IN, $ARGV[0];
my $sum;

while (my $l = <IN>) {
	chomp $l;
	my @a = split /\t/, $l, -1;
	print "$a[0]\t";
	shift @a;
	
	if($ARGV[1] eq "half") {
		my $size = @a;
		for (my $j=0;$j<=$size/2; $j++) {
			shift @a;
		}
	}
	
	$sum  = 0;
	foreach my $i (@a) {
		print "$i\t";
		$sum+=$i;
	}
	print "$sum\n";
}