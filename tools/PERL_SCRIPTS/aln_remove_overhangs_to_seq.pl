BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use ClustalW;
use strict;

my $cl = ClustalW->new;
$cl->setNumCharName(length("031008-E4       "));
$cl->setFile($ARGV[0]);

my $a_ref_aln = $cl->getSeqsWithNames();

my %H = ();
foreach my $s (@$a_ref_aln) {  
  $H{$s->[0]} = $s->[1];
}

my $seq = $H{$ARGV[1]};

my ($gap_to_left)  = $seq =~ /^(\-+)/;
my ($gap_to_right) = $seq =~ /(\-+)$/;


$cl->removeTrailingSequences(length($gap_to_left), 0);
$cl->removeTrailingSequences(length($gap_to_right), 1);

#print length($gap_to_left);
#print "\n";
#print length($gap_to_right);

if (!$ARGV[2]) {
  $cl->delSeq($ARGV[1]);
}

print $cl->getClustalWformat;
