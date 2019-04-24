#!/usr/bin/perl

use lib qw(/home/elemento/PERL_MODULES);

use Table;
use Sets;
use Getopt::Long;

if (scalar(@ARGV) < 2) {
    die "Usage : catcols.pl --file1=ele_bri.txt.evaluation  --file2=ele_bri.txt.evaluation.tfac --cols1=0,1,2,3,4,5,6,7 --cols2=5\n";
}

my $cols      = undef;
my $file1     = undef;
my $file2     = undef;

GetOptions ('cols1=s'      => \$cols1,
	    'cols2=s'      => \$cols2,
	    'file1=s'      => \$file1,
	    'file2=s'      => \$file2);


my @a_cols1 = split /\,/, $cols1;
my @a_cols2 = split /\,/, $cols2;


open IN1, $file1;
open IN2, $file2;


while (($l1 = <IN1>) && ($l2 = <IN2>)) {
    chomp $l1; chomp $l2;
 
    my @a1 = split /\t/, $l1, -1;
    my @a2 = split /\t/, $l2, -1;
    
    my @l = ();

    if ($cols1 eq "all") {
	push @l, @a1;
    } else {
	foreach my $c (@a_cols1) {
	    push @l, $a1[$c];
	}
    }

    if ($cols2 eq "all") {
	push @l, @a2;
    } else {
	foreach my $c (@a_cols2) {
	    push @l, $a2[$c];
	}
    }
    
    print join ("\t", @l) . "\n";
}

close IN1;
close IN2;
