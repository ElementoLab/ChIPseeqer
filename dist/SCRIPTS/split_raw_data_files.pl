#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Sets;
use FileHandle;
use strict;
use Getopt::Long;
use IO::Zlib;

my %FH				= ();
my $num				= 0;
my $format			= "eland"; # can be eland, exteland, mit, bed, sidow, sam 
my $datafolder		= undef;
my $outputfolder	= ".";
my $line			= undef;
my $file			= undef;
my $u				= "U0";
my $verbose         = 1;
my $files           = undef;

#handling lack of arguments
if (@ARGV == 0) {
	die "Lack of arguments. \n Run script as: perl split_raw_data_files.pl --format=STR [ --files=FILES --datafolder=DIR ] --outputfolder=DIR 
	where --format can be eland, exteland, mit, bed, sidow, sam, or export.
	--outputfolder point to a directory where split files will be saved (default = ".", ie current directory)
	--files=FILES, specifies files to process e.g. --files=\"*.gz\"
	
	\n";
}

GetOptions("format=s"	    => \$format,
"verbose=s"      => \$verbose,
"files=s"        => \$files,
"datafolder=s"   => \$datafolder,
"outputfolder=s" => \$outputfolder );

$format = lc($format);

if (!defined($datafolder) & !defined($files)) {
	die "Missing arguments --datafolder or --files. \n Run script as: perl split_raw_data_files.pl --format=STR --datafolder=DIR --outputfolder=DIR \n";
}

if (!defined($outputfolder)) {
	die "Missing argument --outputfolder. \n Run script as: perl split_raw_data_files.pl --format=STR --datafolder=DIR --outputfolder=DIR \n";
}

print STDOUT "Starting splitting files. Keeping uniquely mapping reads. Clonal reads are not filtered out, but will be ignored in the other packages (by default)\n";
print STDOUT "Format is $format\n";
# open data folder
my @contents = ();

if (defined($datafolder)) {
	print STDOUT "Opening directory $datafolder\n";
	opendir (DIRECTORY, $datafolder) or die "can't open directory $datafolder: $!";
	@contents = grep !/^\.\.?$/, readdir DIRECTORY; # skip . and ..
	foreach my $f (@contents) {
		$f = "$datafolder/$f";
	}
} elsif (defined($files)) {
	my $a_ref_f = Sets::getFiles($files);
	@contents = @$a_ref_f;
}

if (! -e $outputfolder) {
	mkdir $outputfolder;
}

# for each file
foreach $file (@contents) {
	next if $file =~ /^reads/;     # skip anything starting with reads
	
	print STDOUT "Opening $file\n";
	
	#open file
	if ($file =~ /\.gz$/) {
		tie *IN, 'IO::Zlib', "$file", "rb";
	} else {
	    open(IN, "$file") or die "can't open file $file: $!";
	} 

	if($format eq "sam") {
	  my $todo = "$ENV{CHIPSEEQERDIR}/split_samfile $file -outdir $outputfolder ";
	  system($todo) == 0 or die "Cannot exec $todo\n";

	} elsif ($format eq "bam") {
	  my $todo = "$ENV{CHIPSEEQERDIR}/split_bamfile $file -outdir $outputfolder ";
	  system($todo) == 0 or die "Cannot exec $todo\n";
	
	} else {
		#for each line
		while ($line = <IN>) {
			
			if($format eq "eland") {
				chomp $line;
				#split line
				my @a = split /[\t]/, $line, -1;
				
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
						$FH{$chr} = new IO::File ">$outputfolder/reads.$chr";
					}
					
					$fh = $FH{$chr};
					
					#write in reads file
					print $fh "$num\t$seq\t$a[2]\t1\t0\t0\t$chr\t$pos\t$str\n"; 
					
					$num++;
				}
			}
			elsif($format eq "exteland") {
				chomp $line;
				
				#split line
				my @a = split /[\t]/, $line, -1;
				if (($a[2] eq '1:0:0') || ($a[2] eq '0:1:0') || ($a[2] eq '0:0:1')) {
					
					#get the sequence
					my $seq = $a[1];
					
					#get the chromosome, chromosome position and strand
					my ($chr, $pos, $str) = $a[3] =~ /(chr.+?)\.fa\:(\d+?)([RF])/;
					
					#create a new reads file
					my $fh = undef;
					
					if (!defined($FH{$chr})) {
						$FH{$chr} = new IO::File ">$outputfolder/reads.$chr";
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
				}
			}
			elsif($format eq "sidow") {
				chomp $line;
				
				#split line
				my @a = split /[\ \:]/, $line, -1;
				
				# get the sequence
				my $seq = $a[0];
				
				#get the chromosome
				my $chr = $a[3];
				
				#get the position
				my $pos = $a[5];
				
				#get the strand
				my $str = $a[6];
				
				#create a new reads file
				my $fh = undef;
				
				if (!defined($FH{$chr})) {
					$FH{$chr} = new IO::File ">$outputfolder/reads.$chr";
				}
				
				$fh = $FH{$chr};
				
				#write in reads file
				print $fh "$num\t$seq\t$u\t1\t0\t0\t$chr\t$pos\t$str\n"; 
				
				$num++;
			}
			elsif(($format eq "bed") || ($format eq "mit")) {
				next if ($line =~ /^track/);
				chomp $line;
				
				#split line
				my @a	= split /[\ \t]/, $line, -1;
				$line	= join("\t", @a);   
				
				#get the chromosome
				my $chr = $a[0];
				my $fh	= undef;
				
				if (!defined($FH{$chr})) {
					$FH{$chr} = new IO::File ">$outputfolder/reads.$chr";
				}
				
				$fh = $FH{$chr};
				
				#write in reads file
				print $fh "$line\n";
				
				$num++;
				
			} 
			elsif($format eq "export") {
				chomp $line;
				
				#split line
				my @a = split /[\t]/, $line, -1;
				
				next if $a[10] !~ /_/;
				
					# get the sequence
					my $seq = $a[8];
					
					#get the chromosome
					my $tmpchr = $a[10];
				
					my @b = split /[_]/, $tmpchr, -1;
					my @c = split /[.]/, $b[2], -1;
					my $chr	= $c[0];
					
					#get the position
					my $pos = $a[12];
					
					#get the strand
					my $str = $a[13];
					
					#create a new reads file
					my $fh = undef;
					
					if (!defined($FH{$chr})) {
						$FH{$chr} = new IO::File ">$outputfolder/reads.$chr";
					}
					
					$fh = $FH{$chr};
					
					#write in reads file
					print $fh "$num\t$seq\t$u\t1\t0\t0\t$chr\t$pos\t$str\n"; 
				
					$num++;
			}
			else { # elsif format
				
				die "Format unknown\n";
			}
			if (($num % 100000 == 0) && ($verbose == 1)) {
				print STDOUT "Wrote $num reads.    \r";
			}
		}  # while (<IN>)
		
	}  # else
	close IN;
}

if (defined($datafolder)) {
	closedir(DIRECTORY);
}

if (($verbose == 1) && ($format ne "sam"))  {
	print STDOUT "Wrote a total of $num reads.\n";
}
