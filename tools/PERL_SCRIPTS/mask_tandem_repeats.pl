use lib "$ENV{HOME}/PERL_MODULES";
use Repeats;
use Sets;
use Fasta;
use strict;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  my $tr = Repeats->new;
  $tr->setSeq($s);
  $tr->setProgram('TRF');
  $tr->process();
  my $a_ref_tandems = $tr->getResults;
  $tr->dispose;
  
  if (@$a_ref_tandems > 0) {
    
  foreach my $r (@$a_ref_tandems) {
    print "$r->[2]\n";
  }

  my $seq_masked = Sets::maskExons($s, $a_ref_tandems, 'X');
  
  print "$s\n";
  print "$seq_masked\n";
  print "\n";
  #<STDIN>;
}
}

