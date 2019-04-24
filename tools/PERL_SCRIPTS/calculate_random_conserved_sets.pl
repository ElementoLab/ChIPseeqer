#
#  the idea is that we run the randomized data generation, calculate a conserved set size for all k-mers and output the 
#

# INPUT : list of k-mers, SEQ1 set, mat.txt, nb of repeats
# OUTPUT : for each k-mer, conserved set
# 

use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


#my $kmerfile = $ARGV[0];
my $seqfile1  = $ARGV[1];
my $seqfile2  = $ARGV[2];
my $seqfile2_random  = $seqfile2; $seqfile2_random =~ s/^.*\/([^\/]+)$/$1/; $seqfile2_random .= ".random";
my $nbrepeats = $ARGV[3];


my $matfile  = $ARGV[4];


my $nbgenes   = `perl ~/PERL_MODULES/SCRIPTS/fasta_sequence_names.pl $seqfile1 | wc -l`; $nbgenes =~ s/[\n\ ]//g;

print "got $nbgenes\n";

if (!defined($matfile)) {
    system("perl ~/PERL_MODULES/SCRIPTS/generate_twospecies_random_datasets.pl --fasta1=$seqfile1 --fasta2=$seqfile2 --outmat=mat.txt --outdir=.");
    $matfile = "mat.txt";
}


my $todo = "";

for (my $i=1; $i<=$nbrepeats; $i++) {

    system("mkdir -p $i/CONSERVED_SETS/");

    
    $todo = "perl ~/PERL_MODULES/SCRIPTS/generate_twospecies_random_datasets.pl --fasta1=$seqfile1 --fasta2=$seqfile2 --freqfile=$matfile --outdir=$i";
    print "$todo\n"; system($todo);

    
    die "$i/$seqfile2_random was not found in $i/\n" if (! -e "$i/$seqfile2_random");

    foreach my $k (@$a_ref) {
	
	$todo = "/home/olly/PROGRAMS/FASTCOMPARE/recompare -re $k->[0] -fasta1 $seqfile1 -fasta2 $i/$seqfile2_random -nbgenes $nbgenes -out out.txt.eval -twostrand 0";
	system($todo);
	
	system("mv out.txt.eval $i/CONSERVED_SETS/$k->[0].txt");
	
    }

    system("rm $i/PSE_D_500.seq.random");

}

