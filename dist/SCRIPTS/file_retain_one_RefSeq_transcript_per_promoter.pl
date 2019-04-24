#!/usr/bin/perl

use lib "$ENV{HOME}/PERL_MODULES";

use Getopt::Long;

use Sets;

my $refseq	= undef;
my $orf		= undef;
my $nm		= undef;
my $chr		= undef;
my $strand	= undef;
my $orfidx		= undef;
my $nmidx		= undef;
my $chridx		= undef;
my $strandidx	= undef;


my %ORFPROM = ();
my %ORFNM   = ();

# handling lack of arguments
if (@ARGV == 0) {
	die "Usage: file_retain_one_RefSeq_transcript_per_promoter.pl --file=FILE --NMidx=INT --ORFidx=INT --chridx=INT --strandidx=INT \n";
}

# processing command line options
GetOptions("file=s"			=> \$refseq,
"NMidx=s"		=> \$nmidx,
"ORFidx=s"		=> \$orfidx,
"chridx=s"		=> \$chridx,
"strandidx=s"	=> \$strandidx);

open IN, $refseq;
my %NM = ();
while (my $l = <IN>) {  
	chomp $l;
	
	my @a = split /\t/, $l, -1;
	
	$orf	= $a[$orfidx];
	$nm		= $a[$nmidx];
	$chr	= $a[$chridx];
	$strand	= $a[$strandidx];
	
	next if ($chr =~ /\_/);
		
	if (defined($ORFNM{$orf}{$nm})) { # we have seen this NM already
		next;
	}
	
	$ORFNM{$orf}{$nm} = 1;
		
	my $tss = undef;
	if ($strand eq "1") {
		$tss = $a[2];
	} else {
		$tss = $a[3];
	}
	
	push @{ $ORFPROM{$orf}{$tss} }, $nm;
	$NM{$nm} = \@a;
	
}
close IN;

foreach my $g (keys(%ORFPROM)) {
		
	foreach my $tss (keys(%{$ORFPROM{$g}})) {
		
		my $small = Sets::get_smallest_NM($ORFPROM{$g}{$tss});
		#print "$g\t$tss\t$small\n";
		#print "$g\t$tss\t$small\t" . join("\t", @{$ORFPROM{$g}{$tss}}) . "\n";		
		print join("\t", @{$NM{$small}}) . "\n";

	}	
}
