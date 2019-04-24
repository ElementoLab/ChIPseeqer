BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use ClustalW;
use strict;

my $cl = ClustalW->new;
$cl->setNumCharName(26);
$cl->setFile($ARGV[0]);

my $a_ref_aln = $cl->getSeqsWithNames();

my @a_n = ();
my @a_s = ();

foreach my $s (@$a_ref_aln) {
  if ($s->[0] =~ /$ARGV[1]/) {
    push @a_n, $s->[0];
    push @a_s, $s->[1];
  }
}

my $cl_out = ClustalW->new;
$cl_out->setNumCharName(26);

$cl_out->setSequences(\@a_n, \@a_s);

print $cl_out->getClustalWformat;

