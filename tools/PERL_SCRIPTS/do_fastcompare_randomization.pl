#
#  the idea is that we run the randomized data generation, calculate a conserved set size for all k-mers and output the 
#

# INPUT : list of k-mers, SEQ1 set, mat.txt, nb of repeats
# OUTPUT : for each k-mer, conserved set
# 

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


my $thehome = $home;
use Table;
use strict;


if (scalar(@ARGV) == 0) {
  
  die "script resfile kmerfile seqfile1 seqfile2 nbrepeats matfile\n";
  
}

my $resfile   = $ARGV[0];
my $kmerfile  = $ARGV[1];
my $seqfile1  = $ARGV[2];
my $seqfile2  = $ARGV[3];
my $seqfile1_random  = $seqfile1; $seqfile1_random =~ s/^.*\/([^\/]+)$/$1/; $seqfile1_random .= ".random";
my $seqfile2_random  = $seqfile2; $seqfile2_random =~ s/^.*\/([^\/]+)$/$1/; $seqfile2_random .= ".random";
my $nbrepeats = $ARGV[4];


my $matfile   = $ARGV[5];
my $out_seq   = "OUT_SEQ";

mkdir $out_seq if (! -e $out_seq);


my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref_res = $ta->getIndex(0);



my $nbgenes = Sets::getNbGenesInFastaFile($seqfile1);

#my $nbgenes   = `perl ~/PERL_MODULES/SCRIPTS/fasta_sequence_names.pl $seqfile1 | wc -l`; $nbgenes =~ s/[\n\ ]//g;



print "got $nbgenes\n";

die "please install ~/PERL_MODULES/SCRIPTS/generate_twospecies_random_datasets_nobioperl.pl\n" if (! -e "$thehome/PERL_MODULES/SCRIPTS/generate_twospecies_random_datasets_nobioperl.pl");

die "please install ~/PROGRAMS/FASTCOMPARE/fastcompare_latest" if (! -e "$thehome/PROGRAMS/FASTCOMPARE/fastcompare_latest");



if (!defined($matfile)) {
  print "matrix file not defined, generate it ..\n";
    system("perl ~/PERL_MODULES/SCRIPTS/generate_twospecies_random_datasets_nobioperl.pl --fasta1=$seqfile1 --fasta2=$seqfile2 --outmat=mat.txt --outdir=.");
    $matfile = "mat.txt";
}


my $todo   = "";
my %SCORES = ();

for (my $i=1; $i<=$nbrepeats; $i++) {
    
    $todo = "perl ~/PERL_MODULES/SCRIPTS/generate_twospecies_random_datasets_nobioperl.pl --fasta1=$seqfile1 --fasta2=$seqfile2 --freqfile=$matfile --outdir=$out_seq --denovo=1";
    print "$todo\n";
    #system($todo);
    print "done\n";

    die "$out_seq/$seqfile1_random was not found in $out_seq/\n" if (! -e "$out_seq/$seqfile1_random");    
    die "$out_seq/$seqfile2_random was not found in $out_seq/\n" if (! -e "$out_seq/$seqfile2_random");

    $todo = "~/PROGRAMS/FASTCOMPARE/fastcompare_latest -kmers $kmerfile -k 7 -nbgenes $nbgenes -fasta1 $out_seq/$seqfile1_random -fasta2 $out_seq/$seqfile2_random > $out_seq/$i.txt";
    print "$todo\n";
    #system($todo);
    print "done\n";

    open IN, "$out_seq/$i.txt";
    while (my $l = <IN>) {
      chomp $l;
      my @a = split /\t/, $l, -1;
      push @{ $SCORES { $a[0] } }, $a[4];
    }
    close IN;
}


foreach my $k ( keys(%SCORES) ) {
  
  next if ($a_ref_res->{ $k }->[ 3 ] < 10);

  my $avg = Sets::average( $SCORES{ $k } );
  my $std = Sets::stddev ( $SCORES{ $k } );

  $std += 0.1;
  my $z   = ($a_ref_res->{ $k }->[ 4 ] - $avg) / $std;

  print "$k\t$z\n";
  
}


