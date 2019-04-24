BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use ClustalW;
use Fasta;
use Sets;
use strict;

# load aa sequences
my $fa_aa1 = Fasta->new;
my $fa_aa2 = Fasta->new;
my $fa_nt1 = Fasta->new;
my $fa_nt2 = Fasta->new;
 
$fa_aa1->setFile($ARGV[0]);
$fa_nt1->setFile($ARGV[1]);

$fa_aa2->setFile($ARGV[2]);
$fa_nt2->setFile($ARGV[3]);




while (1) {

  my $a_ref_aa1 = $fa_aa1->nextSeq();
  my $a_ref_nt1 = $fa_nt1->nextSeq();
  
  my $a_ref_aa2 = $fa_aa2->nextSeq();
  my $a_ref_nt2 = $fa_nt2->nextSeq();
  
  my ($n_aa1, $s_aa1) = @$a_ref_aa1;
  my ($n_nt1, $s_nt1) = @$a_ref_nt1;
  
  my ($n_aa2, $s_aa2) = @$a_ref_aa2;
  my ($n_nt2, $s_nt2) = @$a_ref_nt2;

  

  print $s_aa2 . "\n";
  print $s_nt2 . "\n";
  

  my @a_s = ($s_aa1, $s_aa2);
  my @a_n = ("AA1", "AA2");

  my $cl = ClustalW->new;
  $cl->setSequencesToAlign(\@a_n, \@a_s);
  $cl->run();
    
  my $a_ref_aln = $cl->getInputOrderedSeqArray;
  

  my $s_aa1_aln = shift @$a_ref_aln;
  my $s_aa2_aln = shift @$a_ref_aln;
  
  my $s_nt1_aln = Sets::align_nt_on_aa_sequence($s_aa1_aln, $s_nt1);
  print "$s_nt1_aln\n";

  <STDIN>;

  my $s_nt2_aln = Sets::align_nt_on_aa_sequence($s_aa2_aln, $s_nt2);
  print "$s_nt2_aln\n";
  
  <STDIN>;

}
