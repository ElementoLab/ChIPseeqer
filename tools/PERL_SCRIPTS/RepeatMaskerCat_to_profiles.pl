BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

#
#
# takes as input the .out files from RepeatMasker, and a list of genes
#  build a _index.txt and a names.txt out of it
#
#

use Sets;
use Table;

my $ta = Table->new;
$ta->setDelim("\ +");
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();

shift @$a_ref;
shift @$a_ref;

my %H_NAMES = ();
my %H_INDEX = ();

foreach my $r (@$a_ref) {

  my $n = $r->[5];
  my $s = $r->[10];
  my $t = $r->[11];

  push @{ $H_INDEX{$n} }, $s if (!Sets::in_array($s, @{ $H_INDEX{$n} }));
  push @{ $H_INDEX{$n} }, $t if (!Sets::in_array($t, @{ $H_INDEX{$n} }));
  
  $H_NAMES{ $s } = $s;
  $H_NAMES{ $t } = $t;  
  
  #
  # print "$n\t$s\t$t\n";
  # print join("\t",@$r) . "\n";
  #
  
  

}



open OUT, ">repeat_names.txt";
foreach my $k (keys(%H_NAMES)) {
  print OUT "$k\t$k\tP\n";
}
print OUT "no-repeats\tno-repeats\tP\n";
close OUT;

$ta->loadFile($ARGV[1]);
my $a_ref_seq = $ta->getColumn(0);

open OUT, ">repeat_index.txt";
foreach my $k (@$a_ref_seq) {
  print OUT "$k\t";
  
  if (defined($H_INDEX{$k})) {
    print OUT join("\t", @{$H_INDEX{$k}});
  } else {
    print OUT "no-repeats";
  }

  print OUT "\n";
}
close OUT;
