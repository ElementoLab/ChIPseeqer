#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Sets;
use FileHandle;
use strict;

#handling lack of arguments
if (@ARGV == 0) {
	die "Args: readfiles\n";
}

my %FH	= ();
my $num = 0;
my $u	= "U0";

#for each file
foreach my $f (@ARGV) {
	
	print STDERR "Opening $f\n";
	
	#open file
	open IN, $f;
	
	#for each line
	while (my $l = <IN>) {
		chomp $l;
		
		#split line
		my @a = split /[\t]/, $l, -1;
		
		if (($a[2] eq '1:0:0') || ($a[2] eq '0:1:0') || ($a[2] eq '0:0:1')) {
			
			#get the sequence
			my $seq = $a[1];
			
			#get the chromosome, chromosome position and strand
			my ($chr, $pos, $str) = $a[3] =~ /(chr.+?)\.fa\:(\d+?)([RF])/;
			
			#create a new reads file
			my $fh = undef;
			
			if (!defined($FH{$chr})) {
				$FH{$chr} = new IO::File ">reads.$chr";
			}
			
			$fh = $FH{$chr};
			
			if($a[2] eq '1:0:0') {
				$u="U0";
			}
			elsif($a[2] eq '0:1:0') {
				$u="U1";
			}
			elsif($a[2] eq '0:0:1') {
				$u="U2";
			}
			
			#write in reads file
			print $fh "$num\t$seq\t$u\t1\t0\t0\t$chr\t$pos\t$str\n"; 
			
			$num++;
			
			#print progress
			if (($num % 1000) == 0) {
				print "Wrote $num reads             \r";
			}
		}
	}
	print "\n";
	close IN;
}