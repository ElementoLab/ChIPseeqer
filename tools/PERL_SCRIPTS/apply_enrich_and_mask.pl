BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $todo = "cp $ARGV[1] tmp.fa.1";
system($todo);

for (my $i=1; $i<10; $i++) {
  my $todo = "/home/elemento/PROGRAMS/CLUSTERS/enrich -clusterfile $ARGV[0] -kmerfile /home/elemento/DATA/KMERS/7mers.txt -k 7 -fastafile tmp.fa.$i > enrich.OUT.$i";
  system($todo);

  system("head -1 enrich.OUT.$i");

  my $ta = Table->new;
  $ta->loadFile("enrich.OUT.$i");
  my $a_ref = $ta->getArray();
  my $pat = $a_ref->[0]->[0];

  my $j = $i + 1;
  $todo = "perl /home/elemento/PERL_MODULES/SCRIPTS/mi_mask.pl tmp.fa.$i \"$pat\" > tmp.fa.$j";
  system($todo);
  unlink "tmp.fa.$i";
}
