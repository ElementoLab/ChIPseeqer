#
# actual.txt actual.fa rand.txt rand.fa out.txt out.fa
#

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
use Table;
use Fasta;

open OUTT, ">$ARGV[4]" or die "Cannot open $ARGV[4].\n";
open OUTF, ">$ARGV[5]" or die "Cannot open $ARGV[5].\n";

print OUTT "GENE\tEXP\n";


my $ta = Table->new;

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
foreach my $r (@$a_ref) {
  print OUTT "$r->[0]\t1\n";
}

$ta->loadFile($ARGV[2]);
my $a_ref = $ta->getArray();
foreach my $r (@$a_ref) {
  print OUTT "$r->[0]\t0\n";
}

my $fa = Fasta->new;
$fa->setFile($ARGV[1]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  print OUTF ">$n\n$s\n\n";
  
}

$fa = Fasta->new;

$fa->setFile($ARGV[3]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  print OUTF ">$n\n$s\n\n";
  
}


close OUTF;
close OUTT;
