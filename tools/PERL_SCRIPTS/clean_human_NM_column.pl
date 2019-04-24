BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use strict;

my $genefile = "/home/elemento/PROGRAMS/FIRE/FIRE_DATA/HUMAN/SEQUENCES/human_u_1000_0.fa.genenames";

my $h_ref_genes = Sets::getIndex($genefile);


while (my $l = <STDIN>) {
  chomp $l;
  my @a = split /\ /, $l, -1;

  my $gene = undef;

  foreach my $r (@a) {
    $r =~ s/\.\d+$//;

    #print "'$r'\t$h_ref_genes->{$r}\n";

    if (defined($h_ref_genes->{$r})) {
      $gene = $r;
      last;
    }
    
  }

  if (defined($gene)) {
    print "$gene\n";
  }

}
