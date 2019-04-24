BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use Fasta;
use strict;

my $dir = $ARGV[2];
mkdir $dir if (! -e $dir);

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndexKV(0,1);
my $a_ref_indices = Sets::removeDuplicates($ta->getColumn(1));

foreach my $r (@$a_ref_indices) {

  my $fa = Fasta->new;
  
  $fa->setFile($ARGV[1]);

  my $f = "$dir/$ARGV[0].c$r.seq";

  open OUT, ">$f" or die "cannot open seq file $f\n";
  while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    if ($h_ref->{$n} == $r) {
      print OUT ">$n\n$s\n\n";
    }
    
  } 
  close OUT;


}

