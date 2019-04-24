#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Getopt::Long;

# variables to store the arguments values

my $nm2orf				= undef;
my $genome              = "hg18";		# could be hg18, mm9, dm3, sacser
my $matrix				= undef;
my $mirnafile			= undef;

my %MATRIX				= ();
my @COLS				= ();
my @ROWS				= ();
my %HROWS				= ();
my %NMs					= ();

# handling lack of arguments
if (@ARGV == 0) {
	die "Usage: add_mirna_to_matrix.pl --matrix=FILE --mirnafile=FILE \n";
}

# handling given options
GetOptions("genome=s"			=> \$genome,
"matrix=s"		=> \$matrix,
"mirnafile=s"	=> \$mirnafile);

if (!defined($matrix)) {
	die("Must provide --matrix=FILE\n");
}

if (!defined($mirnafile)) {
	die("Must provide --mirnafile=FILE\n");
}

open MAT, $matrix or die "cannot open $matrix\n";

#open matrix and skip header line
my $l = <MAT>; chomp $l;
my @b = split /\t/, $l, -1;

#store column names
my $n = shift @b;  # ignore
push @COLS, @b;

#read each matrix line
while (my $l = <MAT>) {
	
	chomp $l;	
	my @a = split /\t/, $l, -1;
	my $n = shift @a;
	
	#store the NMs in ROWS
	if (!defined($HROWS{$MIR})) {
		push @ROWS, $n;
	}
	
	$HROWS{$n} = 1;
	
	#store each matrix line in hash (NM)->(values)
	push @{	$MATRIX{$n} }, @a;
}
close MAT;

#new file that stores the matrices and labels
open MATFILES, ">$matrix.MIRfiles.txt";

#open mirnafile
open MIR, $mirnafile or die "cannot open $mirnafile\n";

#read each mirnafile line
while (my $l = <MIR>) {
	
	chomp $l;	
	my @d = split /\t/, $l, -1;
	
	#keep the mirna name
	my $mir = shift @d;

	#for each NM in the mirna
	foreach (@d) {
		#check if NM already exists in matrix
		if (exists($MATRIX{$_})) {
			#store the NM in a hash
			$NMs{$_} += 1;
		}
	}
	
	#for my $key ( keys %NMs ) {
	#	print "$key => $NMs{$key}";
	#	print "\n";
	#}
	
	#print matrix again with new column
	open MATNEW, ">$matrix.$mir.txt";
	
	print MATNEW "GENEID";
	foreach my $c (@COLS) {
		print MATNEW "\t$c";
	}
	print MATNEW "\t$mir\n";
		
	foreach my $g (@ROWS) {
		print MATNEW "$g\t";
		
		print MATNEW join("\t", @{$MATRIX{$g}}); 
		
		if(exists($NMs{$g})) {
			print MATNEW "\t1.0";
		}
		else {
			print MATNEW "\t0.0"; 
		}
		print MATNEW "\n";
	}
	
	close MATNEW;
	
	%NMs = ();
	
	print MATFILES "$matrix.$mir.txt\t$mir\n";
	
}
close MIR;

close MATFILES;