#!/usr/bin/perl
use lib qw(/home/olly/PERL_MODULES);
use ScanACE;
use Table;
use GO_func;
use Yeast;
use Sets;
use Library_Motif_RE;
use GO_categories;
use Hypergeom;


my %h1 = (
	  F1 => "/home/olly/PROGRAMS/FASTCOMPARE/WEBSITE/RESULTS/CER_BAY/DATA/CER.seq",
	  F2 => "/home/olly/PROGRAMS/FASTCOMPARE/WEBSITE/RESULTS/CER_BAY/DATA/BAY.seq",
	  N  => 4358 );

my %h2 = (
	  F1 => "/home/olly/PROGRAMS/FASTCOMPARE/WEBSITE/RESULTS/CER_PAR/DATA/CER.seq",
	  F2 => "/home/olly/PROGRAMS/FASTCOMPARE/WEBSITE/RESULTS/CER_PAR/DATA/PAR.seq",
	  N  => 4695 );
	     
my %h3 = (
	  F1 => "/home/olly/DATA/WORMS/ENSEMBL_ORTHOLOGS/ELE_U.seq",
	  F2 => "/home/olly/DATA/WORMS/ENSEMBL_ORTHOLOGS/BRI_U.seq",
	  N  => 10894 );

my %h4 = (
	  F1 => "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_MOU/HUM_U.seq",
	  F2 => "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_MOU/MOU_U.seq",
	  N  => 15983 );

my %h5 = (
	  F1 => "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_RAT/HUM_U.seq",
	  F2 => "/home/olly/DATA/HUMAN/ENSEMBL_ORTHOLOGS/HUM_RAT/RAT_U.seq",
	  N  => 15386 );

my %h6 = (
	  F1 => "/home/olly/DATA/DROSOPHILA/ORTHOLOGS/OLD/MEL_UPSTREAM.seq",
	  F2 => "/home/olly/DATA/DROSOPHILA/ORTHOLOGS/OLD/PSE_UPSTREAM.seq",
	  N  => 11306 );

my %h7 = ( 
	   F1 => "/home/olly/DATA/YEASTS/ORTHOLOGS/PAR_MIK/PAR_U.seq",
	   F2 => "/home/olly/DATA/YEASTS/ORTHOLOGS/PAR_MIK/MIK_U.seq",
	   N => 4167 );


my %DATA = ( CER_BAY => \%h1,
	     CER_PAR => \%h2,
	     ELE_BRI => \%h3,
	     HUM_MOU => \%h4, 
	     HUM_RAT => \%h5,
	     MEL_PSE => \%h6,
	     PAR_MIK => \%h7
	     );


#
# 
#
my $d1 = $DATA{$ARGV[0]};
my $d2 = $DATA{$ARGV[1]};

#
# load the k-mers to analyze
#
my $t = Table->new;
$t->loadFile($ARGV[2]);
#$t->setLimit(400);
my $a_ref_kmers = $t->getArray;


#
# load the mapping table
#
$t->loadFile($ARGV[3]);
my $a_ref_mapping = $t->getArray();
my $n  = scalar(@$a_ref_mapping);


#
#  nbcopies
#
my $nb = 1;
if ($ARGV[4]) {
	$nb = $ARGV[4];

} else {
	$nb =  1;
	
}

#
# create a mapping orf => num
#
my $k = 1;
my %h_sp1 = ();
my %h_sp2 = ();
foreach my $r (@$a_ref_mapping) {
    $h_sp1{$r->[0]} = $k;
    $h_sp2{$r->[1]} = $k;
    $k++;
}




my $k = 1;
foreach my $r (@$a_ref_kmers) {

    my $re = $r->[0];
    
    $re =~ s/N/\./g;


    #
    #  get conserved set in first species
    #

    
    system("/home/olly/PROGRAMS/FASTCOMPARE/recompare -re \"$re\" -fasta1 $d1->{F1} -fasta2 $d1->{F2} -nbgenes $d1->{N} -out out1.txt -nbcopies $nb > /dev/null");
    my $a_ref1 = Sets::readSet("out1.txt"); 

    
    Sets::writeSet($a_ref1, "TMP/$r->[0]_$ARGV[0].txt");

    #
    #  get conserved set in second species
    #
    my $s_todo = "/home/olly/PROGRAMS/FASTCOMPARE/recompare -re \"$re\" -fasta1 $d2->{F1} -fasta2 $d2->{F2} -nbgenes $d2->{N} -out out2.txt -nbcopies $nb > /dev/null";
    #print "$s_todo\n";
    system($s_todo);

    
    

    my $a_ref2 = Sets::readSet("out2.txt"); 

    Sets::writeSet($a_ref2, "TMP/$r->[0]_$ARGV[1].txt");

    
    #
    #  transform the two sets onto sets with mapped gene names
    #
    
    my @a_set1 = ();
    my @a_set2 = ();

    foreach my $g (@$a_ref1) {
	push @a_set1, $h_sp1{$g} if defined($h_sp1{$g});
    }
    
    foreach my $g (@$a_ref2) {
	push @a_set2, $h_sp2{$g} if defined($h_sp2{$g});
    }


    #
    #  calc significance
    # 
    my $a_ref_ovl = Sets::getOverlapSet(\@a_set1, \@a_set2);
    my $o         = scalar(@$a_ref_ovl);
    my $s1        = scalar(@a_set1);
    my $s2        = scalar(@a_set2);
    
    my $p         = Hypergeom::cumhyper($o, $s1, $s2, $n);

    #
    #  Bonferoni 
    #
    my $nbhyp = scalar(@$a_ref_kmers);
    my $pb    = $p * $nbhyp;

    my $stars = ($pb<0.05?"**":"");
    
    print sprintf("$r->[0]\t$s1\t$s2\t$o\t%e$stars\n", $p);

}


