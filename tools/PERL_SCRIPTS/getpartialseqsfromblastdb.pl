#!/usr/bin/perl
#use lib qw(/home/elemento/PERL_MODULES);
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sequence;
use Sets;

if (!$ARGV[0]) {
	die "Please enter a db and < list of genes ..\n";

}

while (my $l = <STDIN>) {
    
    chomp $l;
    my @a = split /\t/, $l, -1;

    my $o_seq = Sequence->new;

    #$o_seq->setVerbose(1);
    $o_seq->setBlastDB($ARGV[0]);
    
    my $s1 = $o_seq->getSequenceFromBlastDB($a[1], $a[2], $a[3]);


    
    


    if ($s1) {

      if ($a[4] == -1) {
	$s1 = Sets::getComplement($s1);
      }

      print ">$a[0]\n$s1\n\n";
      
      
      
    } else {

	#if ($ARGV[1]) {
	    #print ">$l\n"; print "N" x 100; print "\n\n";	
	#}

    }
    
    
}
