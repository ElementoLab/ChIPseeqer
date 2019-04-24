BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Fasta;
use Alignment;
use Sets;
use strict;

my $al = Alignment->new;
$al->setComp('DNA');
$al->setVerbose(0);
# read seq1
my $fa1 = Fasta->new;
$fa1->setFile($ARGV[0]);
my $a_ref1 = $fa1->nextSeq();
my ($n1, $s1) = @$a_ref1;

# read seq2
my $fa2 = Fasta->new;
$fa2->setFile($ARGV[1]);
my $a_ref2 = $fa2->nextSeq();
my ($n2, $s2) = @$a_ref2;

$s1 = uc($s1);
$s2 = uc($s2);


my $a_ref_aln = $al->sw($s1, $s2);


print $a_ref_aln->[0]->[0];
print "\n";

my $r = $a_ref_aln->[0];
print substr($s1, $r->[1], $r->[2]-$r->[1]);
print "\n\n";

print $a_ref_aln->[1]->[0];
print "\n";

my $r = $a_ref_aln->[1];
print substr($s2, $r->[1], $r->[2]-$r->[1]);
print "\n\n";

