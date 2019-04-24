#!/usr/bin/perl

use lib qw(/home/olly/PERL_MODULES);

use Sequence;
use Sets;
use Table;
use strict;
#
#
# INPUT SP1, SP2, ORT, MIN, MAX
#

if (scalar(@ARGV) < 3) {
    die "Usage : getseqsfromlist.pl DB1 DB2 ORT MIN MAX [OUD] [lenfile deflen] [lenfile2 deflen2]\n";
}

my $min = (defined($ARGV[3])?$ARGV[3]:100);
my $max = (defined($ARGV[4])?$ARGV[4]:1000);

my $thedbfile1  = $ARGV[0];
my $thedbfile2  = $ARGV[1];

my $thegenefile = $ARGV[2];

my $lenfile     = $ARGV[6];
my $deflen      = $ARGV[7];

my $lenfile2     = $ARGV[8];
my $deflen2      = $ARGV[9];

if ($lenfile && !$deflen) {
    die "please provide a default length\n";
}
if ($lenfile2 && !$deflen2) {
    die "please provide a default length\n";
}

open OUT1, ">upstream_sp1.fasta";
open OUT2, ">upstream_sp2.fasta";

my %HH = ();
my $h_ref_len = \%HH;
my %HH2 = ();
my $h_ref_len2 = \%HH;


if (defined($lenfile)) {
    my $t = Table->new;
    $t->loadFile($lenfile);
    $h_ref_len = $t->getIndexKV(0,1);
    
} 

if (defined($lenfile2)) {
    my $t = Table->new;
    $t->loadFile($lenfile2);
    $h_ref_len2 = $t->getIndexKV(0,1);
} 

open INCH, $thegenefile;
while (my $l = <INCH>) {
    
    chomp $l;
 
    my @a = split /\t/, $l;

    my $l1 = $a[0];
    my $l2 = (defined($a[1])?$a[1]:$a[0]);	
	
	    
    my $o_seq = Sequence->new;

    #$o_seq->setVerbose(1);
    $o_seq->setBlastPath("/home/olly/PERL_MODULES/PROGRAMS/BLAST");

    $o_seq->setBlastDB($thedbfile1);    
    my $s1 = $o_seq->getSequenceFromBlastDB("$l1", 0, 0);
    
    $o_seq->setBlastDB($thedbfile2);    
    my $s2 = $o_seq->getSequenceFromBlastDB("$l2", 0, 0);

    if ((length($s1) > $min) && (length($s2) > $min)) {

	#
	#  case where a len is defined (or use the default len);
	#
	if (defined($lenfile) && defined($deflen) && defined($max)) {
	    
	    if ($ARGV[5] eq "D") {
		# length to use 
		my $mylen    = (defined($h_ref_len->{$l1})?$h_ref_len->{$l1}:$deflen);
		my $mylen2   = $mylen;
		if (defined($lenfile2) && defined($deflen2)) {
		    $mylen2    = (defined($h_ref_len2->{$l2})?$h_ref_len2->{$l2}:$deflen2);
		}
		if (defined($max) && ($mylen > $max)) {
		    $mylen = $max;
		}
		if (defined($max) && ($mylen2 > $max)) {
		    $mylen2 = $max;
		}
		
		$s1 = substr($s1, 0, $mylen);
		$s2 = substr($s2, 0, $mylen2);
		
	    } else {
		die "Not implemented yet";
	    }

	}

	elsif ($max) {
            if ($ARGV[5] eq "D") {
            $s1 = substr($s1, 0, $max);
            $s2 = substr($s2, 0, $max);
            } else {
	    $s1 = substr($s1, (length($s1)<$max?0:length($s1) - $max), $max);
	    $s2 = substr($s2, (length($s2)<$max?0:length($s2) - $max), $max);
	    }
	}

	elsif ($ARGV[5] eq "O") {
		$s1 = substr($s1, 3);
		$s2 = substr($s2, 3);
		$s1 = substr($s1, 0, length($s1) - 3);
		$s2 = substr($s2, 0, length($s2) - 3);
	}

	#$l1 =~ s/\.1$//g;
	

	print OUT1 ">$l1\n$s1\n\n";
	print OUT2 ">$l1\n$s2\n\n";

    } else {
	print "$l1\t$l2\t" . length($s1) . "\t" . length($s2) . " (min=$min, $max=$max)\n";
    }	
    
}
    
    
close INCH;
