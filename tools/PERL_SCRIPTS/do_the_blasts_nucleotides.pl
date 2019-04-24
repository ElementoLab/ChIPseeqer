

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
# $mb->setVerbose(1);
$mb->setBlastProgram("blastn");
$mb->setEvalueThreshold("1e-10");
$mb->setNbProcessors(2);

my $a_ref_genomes    = Sets::readSet($ARGV[1]);

my $dir              = $ARGV[2];

my $nbgenomes     = scalar(@$a_ref_genomes);
my $s_tmpstore1 = Sets::getTempFile("/tmp/tmp1.seq");

#
#  traverse all the proteins in file 1


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq()) {
    
    my ($n, $s) = @$a_ref;

    #next if (length($s) < 50);

    
    #
    # get gene name 
    #
    my $gn = $n;

    
    #print "$gn";
    
    #
    # get the protein from DB and save it into a temp file
    #
    open SEQ, ">$s_tmpstore1";
    print SEQ ">$n\n$s\n\n";
    close SEQ;    
        
    # use BLAST to align the current sequence against the reference sequence 
    $mb->setQueryDatabase($s_tmpstore1, "T");

    #
    #  go thru all the genomes
    #

    for (my $i=0; $i<$nbgenomes; $i++) {

      next if (! -e "$dir/$a_ref_genomes->[$i]/whole_genome.fna.nsd");
      
	#die "oops genome file $dir/$a_ref_genomes->[$i]/whole_genome.fna.nsd does not exist ..\n" if (! -e "$dir/$a_ref_genomes->[$i]/whole_genome.fna.nsd");
	
	

	# extract the name of the genome
	my ($g) = ($a_ref_genomes->[$i]); 

	#print "Comparing to $g\n";
	
	# set the database
	$mb->setDatabaseDatabase("$dir/$a_ref_genomes->[$i]/whole_genome.fna", "T");

	
	my $a_hits = $mb->blastallUnique;
	
	if (scalar(@$a_hits) > 0) {
	    #print "gene $gn has a hit in $g\n";
	    print "1\t$g";
	} else {
	    print "0\t$g";
	}

      print "\n";
	
    }

    #print "\n";

}
     
	
unlink $s_tmpstore1;
