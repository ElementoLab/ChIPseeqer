#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
use Fasta;

my %H = ();


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  if (!defined($H{$n})) {
    print ">$n\n$s\n";
    $H{$n} = 1;
  }
}

