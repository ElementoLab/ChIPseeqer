#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Fasta;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $ltot = 0;
while ( my $a_ref = $fa->nextSeq ) {
    my ($name, $seq) = @{$a_ref};
 
    my $l = length($seq);
    #$ltot += $l;
    my @a = split //, $seq;

    foreach my $r (@a) {
	if ($r =~ /[ACGTacgt]/) {
      	  $COUNT { $r } ++;
	  $ltot ++;
        }
    }
}

foreach my $k (keys(%COUNT)) {
  print "$k\t" . ($COUNT{$k}/$ltot) . "\n";
}

my $gc = ($COUNT{"C"}+$COUNT{"G"}) / $ltot;

print "G+C\t$gc\n";


