#!/usr/bin/perl

use lib qw(/home/olly/PERL_MODULES /home/olly/PERL_MODULES/Hypergeom /home/olly/PERL_MODULES/Hypergeom/blib/lib /home/olly/PERL_MODULES/Hypergeom/blib/arch/auto/Hypergeom);

use ScanACE;
use Table;
use GO_func;
use Yeast;
use Sets;
use Library_Motif_RE;
use GO_categories;
use DataFiles;
use strict;
use Getopt::Long;
use File::Copy;
use Insitus;


#
#  parameters
#
my $species      = "CER_BAY";
my $interactions = undef;
my $kmerfile     = undef;
my $condfile     = undef;
my $b_func       = 0;
my $b_insitus       = 0;

my $b_dist       = 0;
my $b_orie       = 0;
my $b_cond       = 0;
my $b_tfac       = 0;
my $b_mots       = 0;
my $b_chip       = 0;
my $verbose      = 0;
my $limit        = undef;
my $after        = 0;
my $createsetsdir = undef;
my $fasta1       = undef;
my $fasta2       = undef;
my $mynbgenes    = undef;
my $bestwins     = undef;
my $typeseq      = 'utr5';
my $pause        = 0;
my $singlestrand    = 0;
my $nbcopies     = undef;
my $myntotal     = undef;

GetOptions ('species=s'        => \$species,
	    'interactions=s'   => \$interactions,
	    'typeseq=s'        => \$typeseq,
	    'kmerfile=s'       => \$kmerfile,    
	    #'condfile=s'       => \$condfile,    
	    'func=i'           => \$b_func,
	    'insitus=i'        => \$b_insitus,
	    'nbcopies=i'       => \$nbcopies,
	    'dist=i'           => \$b_dist,
	    'orie=i'           => \$b_orie,
	    'cond=i'           => \$b_cond,
	    'tfac=i'           => \$b_tfac,
	    'mots=i'           => \$b_mots,
	    'chip=i'           => \$b_chip,
	    'limit=i'          => \$limit,
	    'after=i'          => \$after,
	    'verbose=i'        => \$verbose,
	    'createsetsdir=s'  => \$createsetsdir,
	    'fasta1=s'         => \$fasta1,
	    'fasta2=s'         => \$fasta2,
	    'nbgenes=i'        => \$mynbgenes,
	    'ntotal=i'         => \$myntotal,
	    'pause=i'          => \$pause,
	    'bestwins=s'       => \$bestwins,
	    'singlestrand=i'    => \$singlestrand
	    );


if (!defined($interactions) && !defined($kmerfile)) {
    die "Usage : evaluate_kmers.pl --species=s --typeseq=utr5 --interactions=s --kmerfile=s --condfile=s --func=i --dist=i --orie=i --cond=i --tfac=i --mots=i --chip=i --bestwins=s --verbose=i\n";
    
}

if ($species !~ /CER/) {
    $b_mots = undef;
    $b_chip = undef;
}

my $df      = DataFiles->new;
my $ntotal  = undef;
my $f1      = undef;
my $f2      = undef;
my $nbgenes = undef;
my $allsets = undef;



my $inter   = 0;
if (defined($interactions)) {
    $inter = 1;
}    # 1 to deal with interactions instead

if ($species eq "CER") {
    #$ntotal  = 6000;
    $ntotal  = 5207;
    if ($typeseq eq "utr5") {
	$f1      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_BAY/UPSTREAM/CER_U.seq";
    } else {
	$f1      = "/home/olly/DATA/YEASTS/CEREVISIAE/utr3_cerevisiae_300.fasta";
    }
    
    #$nbgenes = 4358;
}

if ($species eq "ELE") {
    #$ntotal  = 6000;
    $ntotal  = 12832;
    if ($typeseq eq "utr5") {
	$f1      = "";
    } else {
	$f1      = "/home/olly/DATA/WORMS/ENSEMBL_GENES/elegans_downstream_300.fasta";
    }
    
    #$nbgenes = 4358;
}


if ($species eq "MEL") {
    #$ntotal  = 6000;
    $ntotal  = 13525;
    if ($typeseq eq "utr5") {
	$f1      = "";
    } else {
	#$f1      = "/home/olly/DATA/DROSOPHILA/MELANOGASTER/ensembl_droso_300.fasta";
	#$f1      = "/home/olly/DATA/DROSOPHILA/ENSEMBL_GENES/ensembl_droso_3UTR-longest-fbgns.fasta";
	$f1      = "/home/olly/DATA/DROSOPHILA/PSEUDOOBSCURA/MEL_D_500.seq";

    }
    
    #$nbgenes = 4358;
}


if ($species eq "HUM") {
    #$ntotal  = 6000;
    
    if ($typeseq eq "utr5") {
	$f1      = "";
    } else {
	$ntotal  = 21787;
	$f1      = "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/human_downstream_300.fasta";
    }
    
    #$nbgenes = 4358;
}


if ($species eq "CER_BAY") {
    #$ntotal  = 6000;
    $ntotal  = 5207;
    $f1      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_BAY/UPSTREAM/CER_U.seq";
    $f2      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_BAY/UPSTREAM/BAY_U.seq";
    $nbgenes = 4358;
}

if ($species eq "CER_MIK") {
    #$ntotal  = 6000;
    $ntotal  = 5207; # tmp
    $f1      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_MIK/UPSTREAM/CER_U.seq";
    $f2      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_MIK/UPSTREAM/MIK_U.seq";
    $nbgenes = 4308;
}



if ($species eq "CER_KUD") {
    #$ntotal  = 6000;
    $ntotal  = 5207;
    $f1      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_BAY/UPSTREAM/CER_U.seq";
    $f2      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_KUD/UPSTREAM/KUD_U.seq";
    $nbgenes = 4358;
}




if ($species eq "CER_PAR") {
    #$ntotal   = 6000;
    $ntotal   = 5313;
    $f1       = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_PAR/UPSTREAM/CER_U.seq";
    $f2       = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_PAR/UPSTREAM/PAR_U.seq";
    $nbgenes  = 4695;
}

if ($species eq "CER_KLU") {
    #$ntotal  = 6000;
    $ntotal  = 5000; # tmp
    $f1      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_KLU/UPSTREAM/CER_U.seq";
    $f2      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_KLU/UPSTREAM/KLU_U.seq";
    $nbgenes = 2796;
}

if ($species eq "CER_CAS") {
    #$ntotal  = 6000;
    $ntotal  = 5104; # tmp
    $f1      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_CAS/UPSTREAM/CER_U.seq";
    $f2      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_CAS/UPSTREAM/CAS_U.seq";
    $nbgenes = 4113;
}




if ($species eq "ELE_BRI") {
    #$ntotal = 19873;
    $ntotal  = 12832;

    if ($typeseq eq "utr5") {
	$f1      = "/home/olly/DATA/WORMS/ENSEMBL_ORTHOLOGS/ELE_U.seq";
	$f2      = "/home/olly/DATA/WORMS/ENSEMBL_ORTHOLOGS/BRI_U.seq";
	$nbgenes = 10894;
    } else {

	#$f1      = "/home/olly/DATA/WORMS/GENOME/ELE_D_500.seq";
	#$f2      = "/home/olly/DATA/WORMS/GENOME/BRI_D_500.seq";
	#$nbgenes = 10897;
    }

    $allsets = "/home/olly/DATA/WORMS/EXPRESSION/SEPARATE/SETS/allsets.txt";
}

if ($species eq "MEL_PSE") {
    #$ntotal = 19873;
    $ntotal  = 13121;
    $ntotal  = 13522;

    if ($typeseq eq "utr5") {
	$f1      = $df->{"FLIES_MEL_UPSTREAM"};
	$f2      = $df->{"FLIES_PSE_UPSTREAM"};
	$nbgenes = 11306;
    } else {
	$f1      = "/home/olly/DATA/DROSOPHILA/PSEUDOOBSCURA/MEL_D_500.seq";
	$f2      = "/home/olly/DATA/DROSOPHILA/PSEUDOOBSCURA/PSE_D_500.seq";
	
	#$f1      = "/home/olly/DATA/DROSOPHILA/ENSEMBL_GENES/FULL_LENGTH_3UTR/MEL_D_FULL_LENGTH_3000MAX.seq";
	
	#$f2      = "/home/olly/DATA/DROSOPHILA/ENSEMBL_GENES/FULL_LENGTH_3UTR/PSE_D_FULL_LENGTH_3000MAX.seq";


	#$ntotal  = 10263;
	$nbgenes = 11265;

	#$nbgenes = 6752;
       
	
    }

}




if ($species eq "HUM_MOU") {
    #
    #  human / mouse 
    #
    $ntotal  = 8706;
    if ($typeseq eq "utr5") {
	$f1      = "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_MOU/HUM_U.seq";
	$f2      = "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_MOU/MOU_U.seq";
	$nbgenes = 15983;
    } else {
	$f1      = "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_MOU/26JUL2004/HUM_D_300.seq";
	$f2      = "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_MOU/26JUL2004/HUM_D_300.seq";
	$nbgenes = 15977;

    }
  

   
}

if ($species eq "HUM_RAT") {
    #
    #  human / rat
    #
    $ntotal   = 8669;
    $f1       = "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_RAT/HUM_U.seq";
    $f2       = "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_RAT/RAT_U.seq";
    $nbgenes  = 15386;
}


#
#  overload the parameters
#
if (defined($fasta1)) {
    $f1 = $fasta1;
}
if (defined($fasta2)) {
    $f2 = $fasta2;
}
if (defined($mynbgenes)) {
    $nbgenes = $mynbgenes;
}
if (defined($myntotal)) {
    $ntotal = $myntotal;
}

#
#   load the kmers or the interactions
#
my $myfile = (defined($kmerfile)?$kmerfile:$interactions);
my $t = Table->new;
$t->setLimit($limit);
$t->loadFile($myfile);
my $a_ref_kmers = $t->getArray;
print "kmers loaded ..\n" if ($verbose == 1); 

#
#   load the functional data
#
my $y   = undef;
my $m   = undef;
my $goc = undef;

if ($species =~ /CER/) {
    $y = Yeast->new;
    $m = GO_func->new;
    $m->setSource("MIPS", "YEAST");
    $m->setTotalNbORFS($ntotal);
    $m->setPvalueThreshold(0.00001);
    #$m->setBonferroni(1);
}


if ($species =~ /ELE/) {
    #$m = GO_func->new;
    #$m->setSource("GO", "WORM");
    #$m->setTotalNbORFS($ntotal);
    #$m->setPvalueThreshold(0.00001);

    $goc = GO_categories->new;
    $goc->setVerbose(0);
    $goc->setID($df->get("USER"), $df->get("PASS"), $df->get("HOST"));
    $goc->setSpecies("WORMS");
    $goc->setTotalNbORFS($ntotal);
    $goc->setPvalueThreshold(1e-5);
    
    
}

#
#   drosophila
#
if ($species =~ /MEL/) {
    $goc = GO_categories->new;
    $goc->setVerbose(0);
    $goc->setID($df->get("USER"), $df->get("PASS"), $df->get("HOST"));
    $goc->setSpecies("DROSOPHILA");
    $goc->setTotalNbORFS($ntotal);
    #$goc->setTotalNbORFS($nbgenes);
    $goc->setPvalueThreshold(1e-5);
    #$goc->setOriginal(1);

}


#
#   human
#
if ($verbose == 1) {
	print "Connected to DB GO ?\n";
    }
if ($species eq "HUM_MOU") {
    $goc = GO_categories->new;
    $goc->setVerbose(0);
    $goc->setID($df->get("USER"), $df->get("PASS"), $df->get("HOST"));
    $goc->setSpecies("HUMAN");
    $goc->setTotalNbORFS($ntotal);
    $goc->setPvalueThreshold(1.0);
    if ($verbose == 1) {
	print "Connected to DB GO\n";
    }
}


if ($species eq "INT_SAV") {
    $goc = GO_categories->new;
    $goc->setVerbose(0);
    $goc->setID($df->get("USER"), $df->get("PASS"), $df->get("HOST"));
    $goc->setSpecies("CIONA");
    $goc->setTotalNbORFS($ntotal);
    $goc->setPvalueThreshold(0.00001);
    
    if ($verbose == 1) {
	print "Connected to DB GO\n";
    }
}

#
#  load the microarray conditions results
#
if ($species =~ /CER/) {
    $condfile = "/home/olly/RESULTS/MOTIFS_CONDITIONS/YEAST/78mers_nbconditions.txt";
}

if ($species =~ /ELE/) {
    $condfile = "/home/olly/RESULTS/MOTIFS_CONDITIONS/WORM/78mers_nbconditions.txt";
}

if ($species =~ /MEL/) {
    $condfile = "/home/olly/RESULTS/MOTIFS_CONDITIONS/DROSOPHILA/78mers_nbconditions.txt";
}

if ($species =~ /HUM/) {
    $condfile = "/home/olly/RESULTS/MOTIFS_CONDITIONS/HUMAN/78mers_nbconditions.txt";
}




my $h_ref_kmers = undef;
if ($b_cond) {
    $t->setLimit(undef);
    $t->loadFile($condfile);
    $h_ref_kmers = $t->getIndex(0);
}

my $MO = undef;
if ($b_tfac) {
    $MO = Library_Motif_RE->new;
    $MO->loadMotifList("/home/olly/PERL_MODULES/transfac.txt");
    $MO->loadMotifOccurences("/home/olly/PERL_MODULES/TEST/Library_Motif_RE_on_all_kmers.txt");
}

if ($verbose == 1) {
    print "Look at kmers ..\n";
}

my $w = undef;
if (defined($bestwins)) {
    if ($verbose == 1) {
	print "Loading best windows ..\n";
    }
    $t->setLimit(undef);
    $t->loadFile($bestwins);
    $w = $t->getIndex(0);

   
}

my $k = 1;
my $nb_kmers = scalar(@$a_ref_kmers);
my $tmpfile = undef;
foreach my $r (@$a_ref_kmers) {

    if ($k <= $after) {
	$k++;
	next;
    }

    if ($verbose == 1) {
	print "Looking at $r->[0]\n";
    }

    my $re = $r->[0];   
    $re =~ s/N/\./g;
    
    my $re1 = $r->[0];   
    $re1 =~ s/N/\./g;

    my $re2 = $r->[1];   
    $re2 =~ s/N/\./g;

    my $a_ref_dist = undef;
    my $a_ref_pos  = undef;  
    my $a_ref_orie = undef;  
    
    my $med        = undef;
    my $std        = undef;
    my $a_ref      = undef;
    my $siz        = undef;
    my $orie_total = 0;
    my $orie_plus  = 0;

    if ($b_chip || $b_func || $b_mots || $b_dist || $b_orie || $createsetsdir) {
	
	my $s_todo = undef;
	
	if ($inter == 1) {
	    $s_todo = $df->get("MOTIFS_INTERACTIONS") . " -kmer1 $re1 -kmer2 $re2 -fasta1 $f1 -fasta2 $f2 -nbgenes $nbgenes -out out.txt.eval -dist dist.txt.eval"; 
	   
	    if ($singlestrand == 1) {
		$s_todo .= " -twostrand 0";
	    }
	    #print $s_todo;
	    $s_todo .= " > /dev/null";
	} else {

	    if ($species =~ /\_/) {
		$s_todo = $df->get("RECOMPARE") . " -re \"$re\" -fasta1 $f1 -fasta2 $f2 -nbgenes $nbgenes -out out.txt.eval -pos pos.txt.eval -ori ori.txt.eval";
		
		
		if (defined($bestwins) && (defined($w->{$re}->[1])) && (defined($w->{$re}->[2] > 0)) ) {
		    my $m2s = ($w->{$re}->[2] - $w->{$re}->[1]) / 2;
		    my $mmu = $w->{$re}->[1] + $m2s;
		    $s_todo .= " -mu $mmu -2s $m2s";
		}
		
		if ($singlestrand == 1) {
		    $s_todo .= " -twostrand 0";
		}
		
		if (defined($nbcopies)) {
		    $s_todo .= " -nbcopies $nbcopies";
		}

		$s_todo .= " > /dev/null";

		
	    } else {
		
		# single species
		#$s_todo = "/home/olly/PROGRAMS/GENREGEXP/genregexp_singlestrand \"$re\" $f1 | awk '{ print \$1}' | sort | uniq > out.txt.eval";
		# get the full output of genregexp
		$tmpfile = Sets::getTempFile("/tmp/genregexp");
		$s_todo = "/home/olly/PROGRAMS/GENREGEXP/genregexp -re \"$re\" -fastafile $f1 > $tmpfile";

		

		
		
	    }
	}
	
	print "$s_todo\n" if ($verbose == 1);
	system("$s_todo");
	
	if ($species =~ /\_/) { 
	    $a_ref      = Sets::wormRenameSet(Sets::readSet("out.txt.eval")); 
	    
	    #next if (scalar(@$a_ref) == 0);
	    
	    if (defined($createsetsdir)) {
		copy("out.txt.eval", "$createsetsdir/$r->[0].txt");
	    }

	    if ($inter == 1) {
		$a_ref_dist = Sets::readSet("dist.txt.eval"); 
	    } else {
		
		if ($b_dist == 1) {
		    $a_ref_pos  = Sets::readSet("pos.txt.eval");
		}
		
		if ($b_orie == 1) {
		    $a_ref_orie = Sets::readSet("ori.txt.eval"); 
		}
		
	    }

	} else {

	    # ok, read 
	    my $tata = Table->new; 
	    $tata->loadFile($tmpfile);
	    $a_ref = Sets::removeDuplicates($tata->getColumn(0));
	    Sets::writeSet($a_ref, "out.txt.eval");
	    next if (scalar(@$a_ref) == 0);

	    $a_ref_pos = $tata->getColumn(1);
	    #Sets::writeSet($a_ref_pos, "pos.txt.eval");
	    
	    $a_ref_orie = $tata->getColumn(2);
	    #Sets::writeSet($a_ref_ori, "ori.txt.eval");
	    
	    
	    unlink $tmpfile;
	}

	
	
    
	
	#
	# calculate median and standard deviation
	#
	if ($inter == 1) {
	    
	    foreach my $r (@$a_ref_dist) {
		$r = abs($r);
	    }

	    $med = Sets::median($a_ref_dist);
	    $std = Sets::stddev($a_ref_dist);
	} else {

	    if ($b_dist == 1) {
		$med = Sets::median($a_ref_pos);
		$std = Sets::stddev($a_ref_pos);
	    }

	    if ($b_orie == 1) {
		#
		#  count the number of +1
		#
		$orie_total = scalar(@$a_ref_orie);
		foreach my $c (@$a_ref_orie) {
		    if ($c == 1) {
			$orie_plus ++;
		    }
		}
		
		
	    }
	}
    
	
	#
	# size of the ORF set
	#
	$siz = scalar(@$a_ref);

    }

    my $s_func  = undef;
    my $s_prob  = undef;
    my $s_ud    = undef;
    my $s_chip  = undef;
    my $s_tfac  = undef;
    my $s_mot   = undef;
    my $s_conds = undef;
        
    if ($b_func) {

	#
	# mips
	# 

	$s_func = "";
	$s_prob = "";
	
	if ($species =~ /CER/) { 
	    $m->setORFset($a_ref);
	    $m->setVerbose(0);
	    my $a_ref_func = $m->getFunctionalContent;
	    if (scalar(@$a_ref_func) > 0) {
		$s_func = lc($a_ref_func->[0]->{TEXT})
;		$s_prob = sprintf("%d,%3.2e", $a_ref_func->[0]->{S2}, $a_ref_func->[0]->{PVALUE});
		
	    } else {
		$s_func = "";
		$s_prob = "";
	    }
	}


	if ($species eq "HUM_MOU") { 
	    #
	    # human GO
	    #
	    $goc->setEnsemblORFset(Sets::readSet("out.txt.eval", 1));
	    my $a_ref_func = $goc->getFunctionalContent;
	    if (scalar(@$a_ref_func) > 0) {
		$s_func = lc($a_ref_func->[0]->{TEXT});
		$s_prob = sprintf("%3.2e", $a_ref_func->[0]->{PVALUE});
	    } else {
		$s_func = "";
		$s_prob = "";
	    }

	    Sets::writeSet($goc->getORFset(), "out.txt.eval.sp");
	}

	if ($species =~ /ELE/) { 
	    #
	    # drosophila GO
	    #
	    $goc->setORFset(Sets::readSet("out.txt.eval", 1));
	    my $a_ref_func = $goc->getFunctionalContent;
	    if (scalar(@$a_ref_func) > 0) {
		$s_func = lc($a_ref_func->[0]->{TEXT});
		$s_prob = sprintf("%3.2e", $a_ref_func->[0]->{PVALUE});
	    } else {
		$s_func = "";
		$s_prob = "";
	    }
	}


	if (($species =~ /MEL/) || ($species =~ /INT/)) { 
	    #
	    # drosophila GO
	    #
	    $goc->setORFset(Sets::readSet("out.txt.eval", 1));
	    $goc->setMaxCategory(1000);
	    my $a_ref_func = $goc->getFunctionalContent;
	    if (scalar(@$a_ref_func) > 0) {
		$s_func = lc($a_ref_func->[0]->{TEXT});
		#$s_prob = sprintf("%3.2e", $a_ref_func->[0]->{PVALUE});
		$s_prob = sprintf("%d,%d,%d,%d,%3.2e", $a_ref_func->[0]->{S1}, $a_ref_func->[0]->{S2}, $a_ref_func->[0]->{OVERLAP}, $a_ref_func->[0]->{N}, $a_ref_func->[0]->{PVALUE});

	    } else {
		$s_func = "";
		$s_prob = "";
	    }
	}
	
	

    }

    my $s_insitus = "";

    if ($species =~ /MEL/) {
	 my $in = Insitus->new;
	 $in->setORFset(Sets::readSet("out.txt.eval", 1));
	 $in->setBonferroniCorrection(0);
	 
	 $in->setPvalueThreshold(1e-5);
	 $in->calculateEnrichments();
	 my $a_ref_stages_tissues = $in->getStageTissueEnrichments();
	 if (scalar(@$a_ref_stages_tissues)) {
	     my $r1 = shift @$a_ref_stages_tissues;
	     $s_insitus =  sprintf("%s(%3.2e)", $r1->{STAGE_TISSUE}, $r1->{PVALUE}); 

	      
	 }	

	 #print scalar(@$a_ref_stages_tissues); print "\n";
     } 
    
    


    
    #
    # expression conditions
    # 
    if ($b_cond) {

	if (length($r->[0]) > 8) {
	
	    # reduce to 8-mer,
	    my $sg1  = substr($r->[0], 0, length($r->[0]) - 1);
	    my $sg2  = substr($r->[0], 1, length($r->[0]) - 1);
	    my $sg1c = Sets::getComplement($sg1);
	    my $sg2c = Sets::getComplement($sg2);
	    my @a_tmp        = ($h_ref_kmers->{$sg1},  $h_ref_kmers->{$sg2}, 
				$h_ref_kmers->{$sg1c}, $h_ref_kmers->{$sg2c});
	    my @a_tmp_sorted = sort { $a->[1] <=> $b->[1] } @a_tmp;
	    my $a_ref_tmp    = pop @a_tmp_sorted;
	    $s_conds         = $a_ref_tmp->[1];
	    $s_ud            = "$a_ref_tmp->[2]/$a_ref_tmp->[3]";
	    

	} else {
	    
	    if ($h_ref_kmers->{$r->[0]}) {
		$s_conds = $h_ref_kmers->{$r->[0]}->[1];
		$s_ud    = $h_ref_kmers->{$r->[0]}->[2] . "/" .  $h_ref_kmers->{$r->[0]}->[3];
	    } elsif ($h_ref_kmers->{Sets::getComplement($r->[0])}) {
		$s_conds = $h_ref_kmers->{Sets::getComplement($r->[0])}->[1];
		$s_ud    = $h_ref_kmers->{Sets::getComplement($r->[0])}->[2] . "/" .  $h_ref_kmers->{Sets::getComplement($r->[0])}->[3];
	    }
	}
    }
	
    
    #
    # Motifs
    # 
    if ($b_mots) {
	
	my $a_ref_mot = $y->getMotifEnrichmentAndCompareACEscoreForKmer($a_ref, $ntotal, 0.00001, $r->[0], 0.65);
	if (scalar(@$a_ref_mot) > 0) {
	    $s_mot = sprintf("$a_ref_mot->[0]->[0](%3.2e,%3.2f)", $a_ref_mot->[0]->[1], $a_ref_mot->[0]->[2]);
	} else {
	    $s_mot = "";
	}
    }

    #
    # CHIP
    # 
    if ($b_chip) {
	if ($verbose == 1) {
	    print "getting ChIP enrichment ..\n";
	    $y->setVerbose(1);
	}
	my $a_ref_chip = $y->getCHIPEnrichment($a_ref, $ntotal, 0.00001);
	if (scalar(@$a_ref_chip) > 0) {
	    $s_chip = $a_ref_chip->[0]->{REGULATOR} . sprintf("(%3.2e)", $a_ref_chip->[0]->{P});
	} else {
	    $s_chip = "";
	}
    }

    
    #
    # TRANSFAC
    #
    if ($b_tfac) {
	my $txt = "";
	my $a_ref_tfac =  $MO->matchToAll($r->[0]);
	foreach my $t (@$a_ref_tfac) {
	    #print "$t->[0]\t$t->[1]\t$t->[2]\n";
	    my ($name) = $t->[0] =~ /M\d{5}\_(.+?)\.txt/;
	    if ($t->[2] <= 50) { 
		$txt .= "$name($t->[2])/";
	    }
	}
	chop $txt;
	$s_tfac = $txt;
    }

    
    #
    #  
    #
    my $pbino;
    if ($b_orie) {
	
	# test whether the $orie_plus / $orie_total is larger than 0.5
	$pbino = Hypergeom::cumbino($orie_plus, $orie_total, 0.5);

    }
    
    #
    # display the results
    #    
    if ($inter == 1) {
	print sprintf("$k\t%10s\t%10s\t$siz\t", $r->[0], $r->[1]);
	print sprintf("%3.1f", $r->[2]);
    } else {
	print sprintf("$k\t%10s\t", $r->[0]);
	if ($r->[4] =~ /inf/) {
	    print "inf";
	} else {
	    print sprintf("%3.1f", $r->[4]);
	}
    }

    if ($b_dist) {
	#print sprintf("\t$med\t%3.2f", $std);
	print "\t$med";
    }

    if ($b_orie) {
	
	#print sprintf("\t%3.2e%s/%3.2e%s($orie_plus/$orie_total)", 
	#	      $pbino,   ($pbino*$nb_kmers<0.05?"**":""), 
	#	      1-$pbino, ((1-$pbino)*$nb_kmers<0.05?"**":""));
	print "\t$orie_plus\t$orie_total";
    }

    if ($b_cond) {
	$s_ud = "" if ($s_conds eq "");
	print "\t$s_conds $s_ud";
    } 
    
    if ($b_mots) {
	print "\t$s_mot";
    } 

    if ($b_chip) {
	print "\t$s_chip";
    }
    
    if ($b_func) {
	print sprintf("\t%20s", $s_func);
	if ($s_prob) { 
	    print "($s_prob)";
	}
    }

	 if ($b_insitus) {
	print sprintf("\t%20s", $s_insitus);
	#if ($s_prob) { 
	 #   print "($s_prob)";
	#}
    }
	 
    
    if ($b_tfac) {
	my @a_tmp = split /\//, $s_tfac;
	my $a_ref_tmp_nodup = Sets::removeDuplicates(\@a_tmp);
	$s_tfac   = join("/", @$a_ref_tmp_nodup);
	print "\t$s_tfac";
    }
    
    print "\n";

    #system("perl /home/olly/PERL_MODULES/TEST/test_GroupEnrichment.pl out.txt.eval /home/olly/DATA/WORMS/EXPRESSION/SEPARATE/SETS/allsets.txt 10894");
    
    $k++;
    
    <STDIN> if ($pause == 1);

    unlink "out.txt.eval";
    unlink "dist.txt.eval";
    unlink "pos.txt.eval";
}





