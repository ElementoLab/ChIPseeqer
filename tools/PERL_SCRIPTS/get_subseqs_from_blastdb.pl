#!/usr/bin/perl
use lib qw(/home/elemento/PERL_MODULES);

require 'Sequence.pm';

if (!$ARGV[0]) {
	die "Please enter a db and < list of genes ..\n";

}

while (my $l = <STDIN>) {
    
    chomp $l;

    my @a = split /\t/, $l;
    
    my $o_seq = Sequence->new;
    #$o_seq->setVerbose(1);
    $o_seq->setBlastDB($ARGV[0]);
    
    my $s1 = $o_seq->getSequenceFromBlastDB($a[0], $a[3], $a[4]);
    

    if ($s1) {

	print ">$l\n$s1\n\n";
	

	
    } else {

	if ($ARGV[1]) {
	    print ">$l\n"; print "N" x 100; print "\n\n";	
	}

    }
    
    
}
