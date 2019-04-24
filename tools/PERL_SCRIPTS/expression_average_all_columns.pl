#!/usr/bin/perl
use strict;

open IN, $ARGV[0];

#print "ID_REF\tEXP\n";
my $l = <IN>;

print "ID_REF\tEXP\n";

while (my $l = <IN>) {

        
	chomp $l;
	
	my @a = split /\t/, $l, -1;
	my $p = shift @a;

        if ($p eq "ID_REF") {
           print "ID_REF\tlogratio\n";
           next; 
        }

	my $cnt = 0;
	my $sum = 0.0;
	foreach my $r (@a) {
	  if (($r ne "nan") && ($r ne "") && ($r ne "NA")) {
	    $sum += $r;
	    $cnt ++;
	  }
	}
	
	if ($cnt > 0) {
	  my $res = sprintf("%4.3f", $sum / $cnt);
	  print "$p\t$res\n";
	} 
}

close IN;


