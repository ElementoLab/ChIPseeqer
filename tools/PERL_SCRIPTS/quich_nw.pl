BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Alignment;
use Fasta;
use strict;

my $fa1     = Fasta->new;
$fa1->setFile($ARGV[0]);
my $a_ref1  = $fa1->nextSeq(); 
my ($n1, $s1) = @$a_ref1;
$fa1->dispose();


my $fa2     = Fasta->new;
$fa2->setFile($ARGV[1]);
my $a_ref2  = $fa2->nextSeq(); 
my ($n2, $s2) = @$a_ref2;
$fa2->dispose();

my $al = Alignment->new;
$al->setComp('DNA');
#$al->setVerbose(1);

my $a_ref = $al->nw($s1, $s2);

print "$a_ref->[0]\n";
print "$a_ref->[1]\n";

