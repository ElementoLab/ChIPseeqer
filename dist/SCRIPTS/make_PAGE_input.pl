#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Getopt::Long;

# variables to store the arguments values
my $refgene			= undef;
my $geneslist		= undef;
my $geneslistFolder	= undef;
my $outputFolder	= undef;
my $list			= undef;

# handling lack of arguments
if (@ARGV == 0) {
	die "Usage: make_PAGE_input.pl --refgene=FILE --geneslist=FILE [OR --geneslistFolder=FOLDER] \n";
}

# handling given options
GetOptions("refgene=s" => \$refgene,
"geneslist=s"		=> \$geneslist,
"geneslistFolder=s"	=> \$geneslistFolder);

if (!defined($refgene)) {
	die("Must provide --refgene=FILE\n");
}

if (!defined($geneslist)) {
	if(!defined($geneslistFolder)) {
		die("Must provide --genelist=FILE or --geneslistFolder=FOLDER\n");
	}
}

#
# Make RefGene ORF file
#
open(IN, "$refgene") or die "Can't open file $refgene.";

# make new file
open TMP, ">$refgene.ORF_tmp";

while (my $line = <IN>) {
	
	chomp $line;
	
	# split line
	my @a = split /[\t]/, $line, -1;
	
	print TMP "$a[12]\n";
}
close TMP;
close IN;

# set output and list variables
if(defined($geneslist)) {
	$outputFolder	= $geneslist;
	$list			= $geneslist;
}
elsif(defined($geneslistFolder)) {
	$outputFolder	= $geneslistFolder;
	$list			= "$geneslistFolder/allgeneslist.txt"
}

# make ORF file
open ORF, ">$outputFolder\ORF.txt";

if(defined($geneslistFolder)) {
	# open data folder
	my @contents = ();
	
	if (defined($geneslistFolder)) {
		print STDOUT "Opening directory $geneslistFolder\n";
		opendir (DIRECTORY, $geneslistFolder) or die "can't open directory $geneslistFolder: $!";
		@contents = grep !/^\.\.?$/, readdir DIRECTORY; # skip . and ..
		foreach my $f (@contents) {
			print STDOUT "Opening file $f\n";
			$f = "$geneslistFolder/$f";
		}
	} 
	
	# create new file to store all peak files
	open ALLGENESLIST, ">>$list";
	
	# counter for the files
	my $cnt = 1;
	
	# for each file
	foreach $file (@contents) {
		
		next if $file =~ m/allgeneslist/;   # skip anything with allgeneslist
		next if $file =~ m/DS_Store/;		# skip anything with DS_Store
		next if $file =~ m/ORF/;			# skip anything with ORF
		
		print STDOUT "Opening $file\n";
		
		#open file
		open(IN, "$file") or die "can't open file $file: $!";
		
		my %hash = ();
		
		#for each line
		while ($line = <IN>) {
			chomp $line;
			if(!exists($hash{$line})) {
				$hash{$line} +=1;
				print ALLGENESLIST "$line\n";
				print ORF "$line\t$cnt\n"
			}
		}
		$cnt++;
	}
	
	close ALLGENESLIST;
}

#
# Find rest of genes
#
$todo = "$ENV{CHIPSEEQERDIR}/CompareGenes.pl --input1=$refgene.ORF_tmp --input2=$list --output=$outputFolder --keepUniqRec=1";
print "Comparing $list genes list with $refgene.ORF_tmp file...\n";
print("$todo\n");
system($todo) == 0 or die "Cannot exec $todo\n";

if(defined($geneslist)) {
	
	# open genelist file
	open(IN, "$outputFolder\Commons.txt") or die "Can't open file $outputFolder\Commons.txt.";
	
	while (my $line = <IN>) {
		
		chomp $line;
		print ORF "$line\t1\n"
	}
	close IN;
}

# open refgene file
open(IN, "$outputFolder\uunique1.txt") or die "Can't open file $outputFolder\uunique1.txt";

while (my $line = <IN>) {
	
	chomp $line;
	print ORF "$line\t0\n"
}
close IN;


close ORF;

print "File $outputFolder\ORF.txt created.\n";	

print "\nYou are now ready to use iPAGE. Here is a possible command line:\n";

print "page.pl --expfile=$outputFolder\ORF.txt --pathways=human_go_orf\n";