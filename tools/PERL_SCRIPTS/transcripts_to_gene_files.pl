BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use strict;

die "defined U/3/5\n" if (!defined($ARGV[1]));

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my %H = ();
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    my ($g,$t)  = split /\-/, $n;
    my @a_tmp = ($n,$s);
    push @{ $H{ $g} }, \@a_tmp;
}


foreach my $k (keys(%H)) {
  open OUT, ">$ARGV[1]/$k";
  foreach my $r (@{$H{$k}}) {
    print OUT ">$r->[0]\n$r->[1]\n\n";
  }
  close OUT;
}
