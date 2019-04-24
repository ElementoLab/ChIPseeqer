#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Getopt::Long;

# variables to store the arguments values
my $file1			= undef;
my $file2			= undef;
my $output			= undef;

# hashes and hashes-related variables
my %hash1			= ();
my %hash2			= ();

# handling lack of arguments
if (@ARGV == 0) {
	die "Usage: compare_and_join_lists.pl --file1=FILE --file1=FILE \n";
}

# handling given options
GetOptions("file1=s" => \$file1,
"file2=s"		=> \$file2,
"output=s"		=> \$output);

if (!defined($file1)) {
	die("Must provide --file1=FILE\n");
}

if (!defined($file2)) {
	die("Must provide --file2=FILE\n");
}

if (!defined($output)) {
	die("Must provide --output=FILE\n");
}

#
# Open files and put in hashes
#
open(IN, "$file1") or die "Can't open file $file1.";

while (my $line = <IN>) {
	
	chomp $line;
	
	# join line
	my @a = join(" ", split( /[\t]/, $line, -1));
	
	# put in hash
	if(!exists($hash1{$a[0]})) {
		$hash1{$a[0]} +=1;
	}
}
close IN;

open(IN, "$file2") or die "Can't open file $file2.";

while (my $line = <IN>) {
	
	chomp $line;
	
	# join line
	my @a = join(" ", split( /[\t]/, $line, -1));
	
	# put in hash
	if(!exists($hash2{$a[0]})) {
		$hash2{$a[0]} +=1;
	}
}
close IN;

# make new files
#open COMMONS, ">$output\Commons.txt";
#open ALLUNIQUE, ">$output\AllUnique.txt";
#open UNIQUE1, ">$output\uunique1.txt";
#open UNIQUE2, ">$output\uunique2.txt";
open OUT,	">$output\out.txt";

#
# Look for records from file1
#
for my $key ( keys %hash1 ) {
	
	if(exists($hash2{$key})) {
		#print COMMONS "1 \t $key\n";
		print OUT "1\t$key\n";
	}
	else {
		#print UNIQUE1 "0 \t $key\n";
		#print ALLUNIQUE "0 \t $key\n";
		print OUT "0\t$key\n";
	}
}

#
# Look for records from file2
#
for my $key ( keys %hash2 ) {
	
	if(!exists($hash1{$key})) {
		#print UNIQUE2 "0 \t $key\n";
		#print ALLUNIQUE "0 \t $key\n";
		print OUT "0\t$key\n";
	}
}

close OUT;

my $todo = "perl -pi -e 's/ /\t/g' $output\out.txt";
system($todo) == 0 or die "Cannot exec $todo\n";

