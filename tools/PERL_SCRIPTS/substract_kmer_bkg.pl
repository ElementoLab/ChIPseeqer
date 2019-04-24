#!/usr/bin/perl

# prg REAL BKG


# first read bacground
open IN2, $ARGV[1] or die "Cannot open file 2\n";
my %bkg = ();
while (my $l = <IN2>) {
    chomp $l;
    my @a = split /\t/, $l;
    $bkg{ $a[0] } = $a[4];
}
close IN2;


open IN1, $ARGV[0] or die "Cannot open file 1\n";
my @KMERS = ();
while (my $l = <IN1>) {
    #print $l;
    chomp $l;
    my @a = split /\t/, $l;
    $a[5] = ($a[4]=~/inf/?$a[4]:$a[4] - $bkg{ $a[0] });
    push @KMERS, \@a;
    #print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\n";

}
close IN1;

my @SORTED = sort sortInf @KMERS;

foreach my $r (@SORTED) {
    print "$r->[0]\t$r->[1]\t$r->[2]\t$r->[3]\t$r->[5]\n";
}


sub sortInf {
#	my ($a, $b) = @_;
	if ($b->[5] =~ /inf/) {
		return 1;
	} elsif ($a->[5] =~ /inf/) {
		return -1;
	} else {
		return $b->[5] <=> $a->[5];
	}
	
}


