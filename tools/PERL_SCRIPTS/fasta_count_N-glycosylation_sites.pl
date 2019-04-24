#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Fasta;
use strict;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my %COUNT = ();
my $cnt = 0;
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  my $cnts = 0;
  while ($s =~ /(N.[ST])/g) {
    my $d = $1;
    my $p = pos($s);
    #print "$d\t$p\n";
    pos($s) = $p - 3 + 1;
    $cnt ++;
    $cnts ++;
  }
  
  my @a = split /\-/, $n;
  push @{ $COUNT{$a[$#a]} }, $cnts;
  
}

foreach my $g (keys(%COUNT)) {
  my $a =  sprintf("%4.3f", Sets::average($COUNT{$g}));
  print "$g\t$a\tc(" . join(", ", @{$COUNT{$g}}) . ")\n";
}

print "Total: $cnt Glycosylation sites found.\n";
