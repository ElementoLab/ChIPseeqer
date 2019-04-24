#!/usr/bin/perl

#BEGIN{ $home = `echo \$HOME`; chomp $home}
#use lib "$home/PERL_MODULES";
use lib "$ENV{CHIPSEEQERDIR}";

use Sets;
use FileHandle;
use strict;

#handling lack of arguments
if (@ARGV == 0) {
	die "Args: readfiles\n";
}

#my $label = ($ARGV[1] ne ""?$ARGV[1]:"chip");

my %FH = ();
my $num = 0;


#for each file
foreach my $f (@ARGV) {
	
	print STDERR "Opening $f\n";
	
	#open file
	open IN, $f;
	
	#for each line
	while (my $l = <IN>) {
		chomp $l;
		
		#split line
		my @a = split /[\ \:]/, $l, -1;
		
		# get the sequence
		my $seq = $a[0];
		
		#get the chromosome
		my $chr = $a[3];
		
		#get the position
		my $pos = $a[4];
		
		#get the strand
		my $str = $a[5];
		
		#create a new reads file
		my $fh = undef;
		
		if (!defined($FH{$chr})) {
			$FH{$chr} = new IO::File ">reads.$chr";
		}
		
		$fh = $FH{$chr};
		
		#write in reads file
		print $fh "$num\t$seq\tU0\t1\t0\t0\t$chr\t$pos\t$str\n"; 
		
		$num++;
		
		#print progress
		if (($num % 1000) == 0) {
			print "Wrote $num reads             \r";
		}
	}
	print "\n";
	close IN;
}