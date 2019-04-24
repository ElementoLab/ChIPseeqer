#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Getopt::Long;

# variables to store the arguments values

my $nm2orf				= undef;
my $db					= "refSeq";		# could be refSeq, AceView, Ensembl, UCSCGenes
my $genome              = "hg18";		# could be hg18, mm9, dm3, sacser
my $inputfile			= undef;
my $type				= "toNM";		# "toNM" or "toGENE"
my $column				= 0;			# which column to transform
my $outtype				= "replace";    # "replace" a column, or "add" a new column

my %TS2GENE				= ();

# handling lack of arguments
if (@ARGV == 0) {
	die "Usage: transcript_to_gene_transformation.pl --inputfile=FILE --type=\"toNM\" --column=INT \n";
}

# processing command line options
GetOptions("genome=s"			=> \$genome,
"db=s"				=> \$db,
"inputfile=s"		=> \$inputfile,
"type=s"			=> \$type,
"column=s"			=> \$column,
"outtype=s"			=> \$outtype);

$nm2orf	=	"$ENV{CHIPSEEQERDIR}/DATA/$genome/$db.NM2ORF";

open NM2ORF, $nm2orf or die "cannot open $nm2orf\n";
while (my $l = <NM2ORF>) {
	chomp $l;
	my @a = split /\t/, $l, -1;
	
	if($type eq "toNM") {
		$TS2GENE{$a[1]} = $a[0];
	}
	else {
		$TS2GENE{$a[0]} = $a[1];
	}
}
close NM2ORF;

open INPUT, $inputfile or die "cannot open $inputfile\n";
while (my $l = <INPUT>) {
	chomp $l;
	my @b = split /\t/, $l, -1;
	
	if($outtype eq "replace") {
		$l =~ s/$b[$column]/$TS2GENE{$b[$column]}/;
	}
	elsif($outtype eq "add" ) {
		$l .= "\t$TS2GENE{$b[$column]}";
	}
	
	print "$l\n";
}
close INPUT;