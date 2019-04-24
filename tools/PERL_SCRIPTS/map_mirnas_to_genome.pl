use lib qw(/home/olly/PERL_MODULES);

use Fasta;
use Sets;
use MyBlast;
use strict;



my $fa = Fasta->new;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");

$mb->setDatabaseDatabase($ARGV[1]);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("1e-5");
#$mb->setQueryStrand(1);

    
    
#
# go thru all the sequences, align them one by one
#

my $file = Sets::getTempFile("toto");

$fa->setFile($ARGV[0]);

while (my $a_ref_mirna = $fa->nextSeq()) {
  
  my $n = $a_ref_mirna->[0]; $n =~ s/\ .+$//;
  my $s = $a_ref_mirna->[1];
  
  $fa->writeSeq($file, $n, $s);
  $mb->setQueryFile($file);
  
  my $a_ref = $mb->blastallMultiple;
  
  foreach my $hit (@$a_ref) { 
    
    my $hsps = $hit->{"Hit_hsps"};
    
    foreach my $r (@$hsps) {
      
      my $s1 = $r->{"Hsp_align-len"};
      my $e1 = $r->{"Hsp_identity"};
      my $ev = $r->{"Hsp_evalue"};
      my $st = $r->{"Hsp_hit-frame"};
      
      my $hs = $r->{"Hsp_hit-from"};
      my $he = $r->{"Hsp_hit-to"};
      
      
      
      #$n =~ s/\#/\t/;
      
      #next if ($st == -1);
      
      print "$n\t$hit->{Hit_id}\t$hs\t$he\t$s1\t$e1\t$ev\t$st\n";
    }
  }
  
  #<STDIN>;
}

unlink $file;
