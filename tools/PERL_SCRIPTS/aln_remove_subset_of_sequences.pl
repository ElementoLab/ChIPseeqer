BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use ClustalW;
use Sets;
use strict;

my $f = shift @ARGV;

my $cl = ClustalW->new;
$cl->setNumCharName(26);
$cl->setFile($f);

my $a_ref_aln = $cl->getSeqsWithNames();

my @a_n = ();
my @a_s = ();

foreach my $s (@$a_ref_aln) {
  if (!Sets::in_array($s->[0], @ARGV)) {
    push @a_n, $s->[0];
    push @a_s, $s->[1];
  }
}

my $cl_out = ClustalW->new;
$cl_out->setNumCharName(26);

$cl_out->setSequences(\@a_n, \@a_s);

print $cl_out->getClustalWformat;

