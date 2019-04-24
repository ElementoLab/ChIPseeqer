#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use  ClustalW;


my $cl = ClustalW->new;

if (defined($ARGV[1])) {
  $cl->setNumCharName($ARGV[1]);
}

$cl->setFile($ARGV[0]);

my $a_ref_aln = $cl->getSeqsWithNames();

foreach my $s (@$a_ref_aln) {
  print  ">$s->[0]\n$s->[1]\n\n";
}

