#!/usr/bin/perl
use lib qw(/home/elemento/PERL_MODULES);

use Sequence;

use Table;

use Sets;



if (!$ARGV[0]) {
	die "Please enter a list of genes and a DB ..\n";

}

my $o_seq = Sequence->new;
#$o_seq->setVerbose(1);
$o_seq->setBlastDB($ARGV[1]);

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  my $s1 = $o_seq->getSequenceFromBlastDB($r->[1], $r->[2], $r->[3]);
  
  #if (length($s1) >= 100) {
    
    if ($r->[4] == -1) {
      $s1 = Sets::getComplement($s1);
    } 
    
    print ">$r->[0]\n$s1\n\n";
  #}
}

