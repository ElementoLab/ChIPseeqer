BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

$ta->loadFile($ARGV[0]);
my $a_ref_loci = $ta->getColumn(0);
shift @$a_ref_loci;



my $outdir = "$ARGV[1]\_TARGETS";
mkdir $outdir if (! -e $outdir);

foreach my $r (@$a_ref) {
  my $m = shift @$r;
  
  my %H = ();
  foreach my $s (@$r) {
    $H{$s} = 1;
  }

  open OUT, ">$outdir/motif_$m.txt";
  print OUT "GENE\tEXP\n";
  
  foreach my $s (@$a_ref_loci) {
    print OUT "$s\t";
    if (defined($H{$s})) {
      print OUT "1\n";
    } else {
      print OUT "0\n";
    }
  }

  close OUT;

}

