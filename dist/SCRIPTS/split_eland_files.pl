#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Sets;
use FileHandle;
use strict;

if (@ARGV == 0) {
	die "Args: readfiles\n";
}

my %FH = ();
my $num = 0;

print STDERR "Starting splitting files. Keeping U0/U1/U2 reads.\n";

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
		
		if (($a[2] eq 'U0') || ($a[2] eq 'U1') || ($a[2] eq 'U2')) {
			
			#get the sequence, chromosome, chromosome position and strand
			my $seq = $a[1];
			my $chr = $a[6]; $chr =~ s/\.fa//;
			next if ($chr =~ /random/);
			my $pos = $a[7];
			my $str = $a[8];
			
			#create a new reads file
			my $fh = undef;
			
			if (!defined($FH{$chr})) {
				$FH{$chr} = new IO::File ">reads.$chr";
			}
			
			$fh = $FH{$chr};
			
			#write in reads file
			print $fh "$num\t$seq\t$a[2]\t1\t0\t0\t$chr\t$pos\t$str\n"; 
			
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

