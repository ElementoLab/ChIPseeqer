#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Getopt::Long;

# variables to store the arguments values

my $nm2orf				= undef;
my $genome              = "hg18";		# could be hg18, mm9, dm3, sacser
my $inputfile			= undef;

my %TS2MIR				= ();
my @COLS				= ();
my @ROWS				= ();
my %HROWS				= ();

# handling lack of arguments
if (@ARGV == 0) {
	die "Usage: split_mirna_file.pl --inputfile=FILE \n";
}

# handling given options
GetOptions("genome=s"			=> \$genome,
"inputfile=s"		=> \$inputfile);


if (!defined($inputfile)) {
	die("Must provide --inputfile=FILE\n");
}

open INPUT, $inputfile or die "cannot open $inputfile\n";

while (my $l = <INPUT>) {
	chomp $l;
	
	next if ($l =~ m/^\#/);
	
	my @a = split /\t/, $l, -1;
	
	my $NM	= "$a[5]";
	
	$a[1] =~ s/-/_/g;	
	
	my $MIR	= join ("_", $a[0], $a[1]);
	#my $MIR	= "$a[0]";
	
	next if($NM !~ m/^NM/);
	
	if (!defined($HROWS{$MIR})) {
		push @ROWS, $MIR;
	}
	$HROWS{$MIR} = 1;
		
	push @{	$TS2MIR{$MIR} }, $NM;
}
close INPUT;


open OUT, ">$inputfile.MIR.txt";

foreach my $g (@ROWS) {
	print OUT "$g\t";
	print OUT join("\t", @{$TS2MIR{$g}}); 
	print OUT "\n";
}

close OUT;