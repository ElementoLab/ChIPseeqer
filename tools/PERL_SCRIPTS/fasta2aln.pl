#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use ClustalW;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);


my @a_n = ();
my @a_s = ();

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    push @a_n, $n;
    push @a_s, $s;
    
}


my $cl = ClustalW->new;
$cl->setSequences(\@a_n, \@a_s);

print $cl->getClustalWformat;

