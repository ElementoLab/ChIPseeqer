#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use FileHandle;
use strict;


my %MATRIX   = ();
my @COLS     = ();
my $cnt_cols = 0;
my @ROWS     = ();
my %HROWS    = ();

foreach my $f (@ARGV) {
	
	print STDERR "STDERR (not an error): Opening $f\n";
	
	open IN, $f;
	
	# first row is special
	my $l = <IN>; chomp $l;
	my @b = split /\t/, $l, -1;
	
	while (my $l = <IN>) {
		
		chomp $l;
		my @a = split /\t/, $l, -1;
		my $n = shift @a;
		if (!defined($HROWS{$n})) {
			foreach my $c (@COLS) {
				push @{ $MATRIX{$n} }, "";
			}
			push @ROWS, $n;
		}
		$HROWS{$n} = 1;
		
		push @{ $MATRIX{ $n } }, @a;
		
	}
	
	close IN;
	
	# update col num
	my $n = shift @b;  # ignore
	push @COLS, @b;
	#push @COLS, $;
}


print "GENEID";
foreach my $c (@COLS) {
	print "\t$c";
}
print "\n";

foreach my $g (@ROWS) {
	print "$g\t";
	print join("\t", @{$MATRIX{$g}}); 
	print "\n";
}