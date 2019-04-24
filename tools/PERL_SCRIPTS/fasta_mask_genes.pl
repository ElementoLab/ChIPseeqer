#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Sets;
use Table;
if (@ARGV == 0) {
  die "Args: fasta annot\n";
}
my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

my %CHR = ();
foreach my $r (@$a_ref) {
  push @{ $CHR{$r->[1]} }, [ $r->[2], $r->[3] ];
}


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;

    my $s_m = Sets::maskExons($s, $CHR{$n}, 'N');

    $s_m = Fasta::split_large_sequence($s_m, 100);
    
    print ">$n\n$s_m\n";
}

