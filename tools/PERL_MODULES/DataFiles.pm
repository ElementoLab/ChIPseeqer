package DataFiles;


sub new {
    my $self  = {};

    my $home = `echo \$HOME`; 
    chomp $home;
    
    if (($home eq "/") || ($home eq "/root"))  {
	$home = "/home/olly";
    }
    
    $self->{'BLASTDIR'}             = "$home/PERL_MODULES/PROGRAMS/BLAST/bin";
    $self->{'BLASTALL'}             = $self->{'BLASTDIR'} . "/blastall";
    $self->{'FASTACMD'}             = $self->{'BLASTDIR'} . "/fastacmd";
    $self->{'GENREGEXP'}            = "$home/PROGRAMS/GENREGEXP/genregexp";

    $self->{'RECOMPARE'}            = "$home/PROGRAMS/FASTCOMPARE/recompare";
    $self->{'KMERS_LOCATE'}         = "$home/PROGRAMS/KMERS_LOCATE/locate_kmers";
    $self->{'FASTA33'}              = "$home/PERL_MODULES/PROGRAMS/FASTA/fasta33";
    $self->{'TFASTX33'}              = "$home/PERL_MODULES/PROGRAMS/FASTA/tfastx33";
    $self->{'SCANACE'}              = "$home/PERL_MODULES/PROGRAMS/ACE/ScanACE";

    $self->{'MOTIFS_INTERACTIONS'}  = "$home/PROGRAMS/MOTIFS_INTERACTIONS/motifs_interaction";

    $self->{'REPEATMASKER'}         = "$home/PERL_MODULES/PROGRAMS/REPEATMASKER/RepeatMasker/RepeatMasker";
    $self->{'TRF'}                  = "$home/PROGRAMS/SOFT/trf404.mac-leopard.exe";
    $self->{'LALIGN'}               = "$home/PERL_MODULES/PROGRAMS/LALIGN/lalign";
    $self->{'DIALIGN'}              = "$home/PERL_MODULES/PROGRAMS/DIALIGN/dialign2-2";

    #
    #  fastcompare data
    #
    $self->{'WORMS_ELE_UPSTREAM'}   = "$home/DATA/WORMS/ENSEMBL_ORTHOLOGS/ELE_U.seq";
    $self->{'WORMS_ELE_DOWNSTREAM'} = "$home/DATA/WORMS/ENSEMBL_ORTHOLOGS/ELE_D.seq";
    $self->{'FLIES_MEL_UPSTREAM'}   = "$home/DATA/DROSOPHILA/ORTHOLOGS/MEL_U.seq";
    $self->{'FLIES_MEL_DOWNSTREAM'} = "$home/DATA/DROSOPHILA/ORTHOLOGS/MEL_D.seq";
    $self->{'FLIES_PSE_UPSTREAM'}   = "$home/DATA/DROSOPHILA/ORTHOLOGS/PSE_U.seq";
    $self->{'FLIES_PSE_DOWNSTREAM'} = "$home/DATA/DROSOPHILA/ORTHOLOGS/PSE_D.seq";

    
    #
    #  fastcompare results
    #
    $self->{'FASTCOMPARE_CER_BAY_CLUSTERED'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_BAY/UPSTREAM/CLUSTERED/cer_bay.txt.withsubstrings";
    $self->{'FASTCOMPARE_CER_BAY_CLUSTERED'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_PAR/UPSTREAM/CLUSTERED/cer_par.txt.withsubstrings";

    $self->{'FASTCOMPARE_CER_BAY_NONCLUSTERED_7MERS'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_BAY/UPSTREAM/NONCLUSTERED/7mers_all.txt";
    $self->{'FASTCOMPARE_CER_BAY_NONCLUSTERED_8MERS'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_BAY/UPSTREAM/NONCLUSTERED/8mers_all.txt";
    $self->{'FASTCOMPARE_CER_BAY_NONCLUSTERED_9MERS'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_BAY/UPSTREAM/NONCLUSTERED/9mers_all.txt";
    $self->{'FASTCOMPARE_CER_BAY_NONCLUSTERED_10MERS'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_BAY/UPSTREAM/NONCLUSTERED/10mers_all.txt";


    $self->{'FASTCOMPARE_CER_PAR_NONCLUSTERED_7MERS'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_PAR/UPSTREAM/7mers_all.txt";
    $self->{'FASTCOMPARE_CER_PAR_NONCLUSTERED_8MERS'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_PAR/UPSTREAM/8mers_all.txt";
    $self->{'FASTCOMPARE_CER_PAR_NONCLUSTERED_9MERS'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_PAR/UPSTREAM/9mers_all.txt";
    $self->{'FASTCOMPARE_CER_PAR_NONCLUSTERED_10MERS'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_PAR/UPSTREAM/10mers_all.txt";



    #
    #  fastcompare corrected results
    #
    $self->{'FASTCOMPARE_CER_BAY_CLUSTERED_CORRECTED'}            = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_BAY/UPSTREAM/CORRECTED/cer_bay.txt.withsubstrings";
    $self->{'FASTCOMPARE_CER_BAY_NONCLUSTERED_7MERS_CORRECTED'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_BAY/UPSTREAM/CORRECTED/7mers_all.txt";
    $self->{'FASTCOMPARE_CER_BAY_NONCLUSTERED_8MERS_CORRECTED'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_BAY/UPSTREAM/CORRECTED/8mers_all.txt";
    $self->{'FASTCOMPARE_CER_BAY_NONCLUSTERED_9MERS_CORRECTED'}   = "$home/RESULTS/FASTCOMPARE/YEASTS/CER_BAY/UPSTREAM/CORRECTED/9mers_all.txt";


    #
    #  GO
    #
    $self->{'GO_PROCESS_ONTOLOGY'}  = "$home/PERL_MODULES/FUNCTIONS/GO/ONTOLOGY/process.ontology";
    $self->{'GO_FUNCTION_ONTOLOGY'}  = "$home/PERL_MODULES/FUNCTIONS/GO/ONTOLOGY/function.ontology";
    $self->{'GO_COMPONENT_ONTOLOGY'}  = "$home/PERL_MODULES/FUNCTIONS/GO/ONTOLOGY/component.ontology";

    #
    #  known motifs
    #
    $self->{'KNOWN_MOTIFS_LIST'}    = "$home/DATA/YEASTS/KNOWN_MOTIFS/known_motifs.txt";

    #
    #  yeast 1kb upstream sequences
    #
    $self->{'CER_UPSTREAM'}   = "$home/DATA/YEASTS/CEREVISIAE/utr5_sc_1000.fasta";
    $self->{'PAR_UPSTREAM'}   = "$home/DATA/YEASTS/PARADOXUS_MIT/utr5_1000_paradoxus_mit.fasta";
    $self->{'BAY_UPSTREAM'}   = "$home/DATA/YEASTS/BAYANUS_MIT/utr5_1000_bayanus_mit.fasta";
    $self->{'KUD_UPSTREAM'}   = "$home/DATA/YEASTS/KUDRIAVZEVII_WASHU/utr5_1000_kudriavzevii_washu.fasta";
    $self->{'MIK_UPSTREAM'}   = "$home/DATA/YEASTS/MIKATAE_MIT/utr5_1000_mikatae_mit.fasta";
  

    #
    # yeast multiple alignments
    #
    $self->{'YEAST_ALIGNMENTS'} = "$home/DATA/YEASTS/ALIGNMENTS/list.txt";
    

    $self->{'ELE_DOWNSTREAM_300'}   = "$home/DATA/WORMS/ENSEMBL_GENES/elegans_downstream_300.fasta";

    


    #
    #  
    #
    $self->{'YEAST_ORFS_GENES'} = "$home/DATA/YEASTS/CEREVISIAE/cerevisiae_orfs_genes.txt";

    #
    #  KMERS
    #
    $self->{'7MERS'}  = "/home/olly/DATA/KMERS/7mers.txt";
    $self->{'8MERS'}  = "/home/olly/DATA/KMERS/8mers.txt";
    $self->{'9MERS'}  = "/home/olly/DATA/KMERS/9mers.txt";
    $self->{'10MERS'} = "/home/olly/DATA/KMERS/10mers.txt";
    


    #
    #   bi-species stuff
    #

    #
    #   YEAST
    #
    my %h_tmp = (
		 'LEN'     => 1000,	       
		 "FASTA1"  => "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_BAY/UPSTREAM/CER_U.seq",
		 "FASTA2"  => "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_BAY/UPSTREAM/BAY_U.seq",
		 "NBGENES" => 4358,
		 "NBGO"    => 5207
		 );
    $self->{'CER_BAY'} = \%h_tmp;
        
    
    #
    #   FLIES
    #
    my %h_tmp = (
		 'LEN'     => 2000,	       
		 "FASTA1"  => "/home/olly/DATA/DROSOPHILA/ORTHOLOGS/MEL_U.seq",
		 "FASTA2"  => "/home/olly/DATA/DROSOPHILA/ORTHOLOGS/PSE_U.seq",
		 "NBGENES" => 11306,
		 "NBGO"    => 13121
		 );
    $self->{'MEL_PSE'} = \%h_tmp;


    #
    #   WORMS
    #
    my %h_tmp = (
		 'LEN'     => 2000,	       
		 "FASTA1"  => "/home/olly/DATA/WORMS/ENSEMBL_ORTHOLOGS/ELE_U.seq",
		 "FASTA2"  => "/home/olly/DATA/WORMS/ENSEMBL_ORTHOLOGS/BRI_U.seq",
		 "NBGENES" => 10894,
		 "NBGO"    => 12832 
		 );
    $self->{'ELE_BRI'} = \%h_tmp;


    #
    #   HUMAN
    #
    my %h_tmp = (
		 'LEN'     => 2000,	       
		 "FASTA1"  => "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_MOU/HUM_U.seq",
		 "FASTA2"  => "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_MOU/MOU_U.seq",
		 "NBGENES" => 15983,
		 "NBGO"    => 8706
		 );
    $self->{'HUM_MOU'} = \%h_tmp;

    

    #
    # password and stuff
    #
    my $h = `echo \$HOSTNAME`; chomp $h;
    
    #if (($h eq "gen-laptop-elemento.Princeton.EDU") || ($h eq "localhost.localdomain")) {
	$self->{'USER'} = "root";
	$self->{'PASS'} = "";
	$self->{'HOST'} = "localhost";
    #} else {
#	$self->{'USER'} = "olly";
#	$self->{'PASS'} = "batman";
#	$self->{'HOST'} = "128.112.116.214";
    #}
    
    #print $self->{'HOST'};

    bless($self);
    return $self;
}



sub get {
   
    my ($self, $s) = @_;
    
    return $self->{$s};

}


1;
