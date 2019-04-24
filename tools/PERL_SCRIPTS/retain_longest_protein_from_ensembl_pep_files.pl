#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Fasta;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while ( my $a_ref = $fa->nextSeq ) {
    my ($name, $seq) = @{$a_ref};

    my ($gene) = $name =~ /gene\:(.+?)\ /;
    
    if (!defined($SEQ{$gene})) {
      $SEQ{$gene} = $seq;
    } elsif (length($SEQ{$gene}) < length($seq)) {
      $SEQ{$gene} = $seq;
    }
    
    #print "$name\n";
}


foreach my $g (sort(keys(%SEQ))) {
  print ">$g\n$SEQ{$g}\n\n";
}
