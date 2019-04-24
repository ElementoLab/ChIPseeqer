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

#
#  parameters
#
my $species      = "CER_BAY";
my $interactions = undef;
my $kmerfile     = undef;
my $condfile     = undef;
my $b_func       = 1;
my $b_dist       = 1;
my $b_cond       = 0;
my $b_tfac       = 0;
my $b_mots       = 0;
my $b_chip       = 0;
my $verbose      = 0;
my $limit        = undef;
my $after        = 0;
my $dir          = undef;
my $nbcopies     = 1;
GetOptions ('species=s'        => \$species,
	    'interactions=s'   => \$interactions,
	    'kmerfile=s'       => \$kmerfile,    
	    'limit=i'          => \$limit,
	    'after=i'          => \$after,
	    'verbose=i'        => \$verbose,
	    'dir=s'            => \$dir,
	    'nbcopies=s'       => \$nbcopies
	    );

if (!defined($interactions) && !defined($kmerfile)) {
    die "Usage : generate_conserved_sets.pl --species=s --interactions=s --kmerfile=s --dir=s --nbcopies X\n";
    
}


my $df      = DataFiles->new;
my $ntotal  = undef;
my $f1      = undef;
my $f2      = undef;
my $nbgenes = undef;

my $inter   = 0;
if (defined($interactions)) {
    $inter = 1;
}    # 1 to deal with interactions instead


if ($species eq "CER_BAY") {
    #$ntotal  = 6000;
    $ntotal  = 5207;
    $f1      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_BAY/UPSTREAM/CER_U.seq";
    $f2      = "/home/olly/DATA/YEASTS/ORTHOLOGS/CER_BAY/UPSTREAM/BAY_U.seq";
    $nbgenes = 4358;
}




if ($species eq "CER_PAR") {
    #$ntotal   = 6000;
    $ntotal   = 5313;
    $f1       = "../FASTCOMPARE/WEBSITE/RESULTS/CER_PAR/DATA/CER.seq";
    $f2       = "../FASTCOMPARE/WEBSITE/RESULTS/CER_PAR/DATA/PAR.seq";
    $nbgenes  = 4695;
}


if ($species eq "PAR_MIK") {
    #$ntotal   = 6000;
    $ntotal   = 5313;
    $f1       = "/home/olly/DATA/YEASTS/ORTHOLOGS/PAR_MIK/PAR_U.seq";
    $f2       = "/home/olly/DATA/YEASTS/ORTHOLOGS/PAR_MIK/MIK_U.seq";
    $nbgenes  = 4167;
}


if ($species eq "ELE_BRI") {
    #$ntotal = 19873;
    $ntotal  = 12832;
    $f1      = "/home/olly/DATA/WORMS/ENSEMBL_ORTHOLOGS/ELE_U.seq";
    $f2      = "/home/olly/DATA/WORMS/ENSEMBL_ORTHOLOGS/BRI_U.seq";
    $nbgenes = 10894;
}

if ($species eq "MEL_PSE") {
    #$ntotal = 19873;
    $ntotal  = 13121;
    $f1      = $df->{"FLIES_MEL_UPSTREAM"};
    $f2      = $df->{"FLIES_PSE_UPSTREAM"};
    $nbgenes = 11306;
}


if ($species eq "HUM_MOU") {
    #
    #  human / mouse 
    #
    $ntotal  = 8706;
    $f1      = "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_MOU/HUM_U.seq";
    $f2      = "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_MOU/MOU_U.seq";
    $nbgenes = 15983;
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
#   load the kmers or the interactions
#
my $myfile = (defined($kmerfile)?$kmerfile:$interactions);
my $t = Table->new;
$t->setLimit($limit);
$t->loadFile($myfile);
my $a_ref_kmers = $t->getArray;
print "kmers loaded ..\n" if ($verbose == 1); 


my $k = 1;
foreach my $r (@$a_ref_kmers) {

    if ($k <= $after) {
	$k++;
	next;
    }

    my $re  = $r->[0];   
    $re     =~ s/N/\./g;
    
    my $re1 = $r->[0];   
    $re1    =~ s/N/\./g;

    my $re2 = $r->[1];   
    $re2    =~ s/N/\./g;

    my $s_todo = undef;
    
    if ($inter == 1) {
	$s_todo = $df->get("MOTIFS_INTERACTIONS") . " -kmer1 $re1 -kmer2 $re2 -fasta1 $f1 -fasta2 $f2 -nbgenes $nbgenes -maxov 0 -out out.txt -dist dist.txt";
    } else {
	$s_todo = $df->get("RECOMPARE") . " -re \"$re\" -fasta1 $f1 -fasta2 $f2 -nbgenes $nbgenes -out out.txt -pos pos.txt";

	if ($nbcopies > 1) {
	    $s_todo .= " -nbcopies $nbcopies";
	}
    }

    print "$s_todo\n" if ($verbose == 1);
    system("$s_todo > /dev/null");
    
    my $a_ref      = Sets::readSet("out.txt"); 

    my $filename   = undef;
    if ($inter == 1) {
	$filename  = $re1 . "_" . $re2  . ".txt";
    } else {
	$filename  = "$re.txt";
    }
    Sets::writeSet($a_ref, "$dir/$filename");


}

