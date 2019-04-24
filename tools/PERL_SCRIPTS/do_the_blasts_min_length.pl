

#
#  INPUT : one set of aa sequences, one list of genomes
#

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use lib "$home/usr/lib/perl5/site_perl/5.6.1";
use lib "$home/usr/lib/perl5/site_perl/5.8.3";

use MyBlast;
use Fasta;
use Sets;
use Sequence;
use Table;
use DataFiles;
use strict;

my $df = DataFiles->new;

#
#  get a new sequence object to retrieve BLAST seqs
#


my $mb = MyBlast->new;
#$mb->setVerbose(1);
$mb->setBlastProgram("blastp");
$mb->setEvalueThreshold("1e-10");
$mb->setNbProcessors(1);

my $a_ref_genomes    = Sets::readSet($ARGV[1]);

my $nbgenomes     = scalar(@$a_ref_genomes);
my $s_tmpstore1 = Sets::getTempFile("/tmp/tmp1.seq");

# extract the nameof the genome
foreach my $gg (@$a_ref_genomes) {
  my ($g) = $gg =~ /GENOMES\/(.+?)\/genome/;
  print "\t$g";
}
print "\n";

#
#  traverse all the proteins in file 1
#

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
    
    my ($n, $s) = @$a_ref;

    next if (length($s) < 50);

    
    #
    # get gene name 
    #
    my ($gn) = $n =~ /gi\|\d+\|ref\|(.+?)\|/;


    #print "$n\n";

    print "$gn";
    
    #
    # get the protein from DB and save it into a temp file
    #
    $s =~ s/\r//g;
    open SEQ, ">$s_tmpstore1";
    print SEQ ">$n\n$s\n\n";
    close SEQ;
    
        
    # use BLAST to align the current sequence against the reference sequence 
    $mb->setQueryDatabase($s_tmpstore1, "T");


    #
    #  go thru all the genomes
    #

    for (my $i=0; $i<$nbgenomes; $i++) {

	next if (! -e "$a_ref_genomes->[$i].psd");
	
	

	# extract the nameof the genome
	my ($g) = $a_ref_genomes->[$i] =~ /GENOMES\/(.+?)\/genome/;

	#print "Comparing to $g\n";
	
	# set the database
	$mb->setDatabaseDatabase($a_ref_genomes->[$i], "T");

	
	my $a_ref_hsps = $mb->blastallUnique;
	
	if (scalar(@$a_ref_hsps) > 0) {
	  my $query_length = $mb->getQueryLength();
	  my $hit_length   = $mb->getUniqueHitLength();
	  my @a_hsps_q = ();
	  my @a_hsps_d = ();
	  foreach my $h (@$a_ref_hsps) {
	    
	    #print "raw Q: $h->{QFROM} => $h->{QTO}\n";
	    #print "raw D: $h->{DFROM} => $h->{DTO}\n";

	    my @b = ($h->{"QFROM"}, $h->{"QTO"});
	    push @a_hsps_q, \@b;
	    my @c = ($h->{"DFROM"}, $h->{"DTO"});
	    push @a_hsps_d, \@c;
	  }

	  my $a_ref_q = Sets::assemble_overlapping_fragments(\@a_hsps_q);
	  my $a_ref_d = Sets::assemble_overlapping_fragments(\@a_hsps_d);
	  
	  my $frac_identical1 = 0;
	  foreach my $f (@$a_ref_q) {
	    #print "Q: got $f->[0] -> $f->[1]\n";
	    $frac_identical1 += ( $f->[1] - $f->[0] + 1);
	  }
	  
	  my $frac_identical2 = 0;
	  foreach my $f (@$a_ref_d) {
	    #print "D: got $f->[0] -> $f->[1]\n";
	    $frac_identical2 += ( $f->[1] - $f->[0] + 1);
	  }
	  
	  #
	  my $frac_len_seq1    = $frac_identical1 / $query_length;
	  my $frac_len_seq2    = $frac_identical2 / $hit_length;

	  my $frac = 0.5;
	  if (($frac_len_seq1 > $frac) && ($frac_len_seq2 > $frac)) {
	    print "\t1";
	  } else {
	    print "\t0";
	  }
	  
	  
	} else {
	  print "\t0";
	}
	
	
	
      }
  
    print "\n";
    
  }
     
	
unlink $s_tmpstore1;
