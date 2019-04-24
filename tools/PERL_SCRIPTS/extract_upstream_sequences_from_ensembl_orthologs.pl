#!/usr/bin/perl

use strict;
use Bio::Seq;
use lib qw(/home/olly/PERL_MODULES);
use Sequence;
use Sets;
use DataFiles;


my $thedbfile1 = $ARGV[2];
my $thedbfile2 = $ARGV[3];


my $df = DataFiles->new;

if (scalar(@ARGV) == 0) {
    print "Usage : THISPRG predictions1 predictions2 DB1 DB2 min max\n";
    exit;
}


my $min        = undef;
if (defined($ARGV[4])) { 
    $min        = $ARGV[4];
} else {
    $min        = 100;
}

my $max        = undef;
if (defined($ARGV[5])) { 
    $max        = $ARGV[5];
} else {
    $max        = 10000000000;
}

my $h_ref1 = Sets::readEnsemblOrthologs($ARGV[0]);

my $h_ref2 = Sets::readEnsemblOrthologs($ARGV[1]);

my $a_ref  = Sets::getReciprocalOrthologs($h_ref1, $h_ref2);


#
#   extract  the data
#

open OUT1, ">upstream_sp1.fasta" or die "Cannot open upstream_sp1.fasta for writing ..\n";
open OUT2, ">upstream_sp2.fasta" or die "Cannot open upstream_sp2.fasta for writing ..\n";


my $cnt = 0;
foreach my $r (@$a_ref) {
    
    if ($ARGV[0] =~ /ele/) { 
      $r->[0] .= ".1";
    }
    
    print "$r->[0]\t$r->[1]\n";
    
    my $o_seq = Sequence->new;
    
    $o_seq->setBlastPath($df->get("BLASTDIR"));
    
    $o_seq->setBlastDB($thedbfile1);    
    my $s1 = $o_seq->getSequenceFromBlastDB($r->[0], 0, 0);
    
    $o_seq->setBlastDB($thedbfile2);    
    my $s2 = $o_seq->getSequenceFromBlastDB($r->[1], 0, 0);

    
    #print "s1=$s1\ts2=$s2\n";

    if ((length($s1) > $min) && (length($s2) > $min)) {

	$s1 = substr($s1, 0, $max);
	$s2 = substr($s2, 0, $max);
		
	print OUT1 ">$r->[0]\n$s1\n\n";
	print OUT2 ">$r->[0]\n$s2\n\n";

    }	
    
    $cnt++;

    #last if ($cnt == 10);
}
    


close OUT1;
close OUT2;







