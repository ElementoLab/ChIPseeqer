use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $seqfile1 = $ARGV[1];
my $seqfile2 = $ARGV[2];
my $dir      = $ARGV[4];
my $nbgenes  = $ARGV[3]; 


foreach my $k (@$a_ref) {
    
    $todo = "/home/olly/PROGRAMS/FASTCOMPARE/recompare -re $k->[0] -fasta1 $seqfile1 -fasta2 $seqfile2 -nbgenes $nbgenes -out out.txt.eval -twostrand 0";
    system($todo);
    
    system("mv out.txt.eval $dir/$k->[0].txt");

    #print "Conserved sets generated for $k->[0]\n";
}

#system("rm $i/PSE_D_500.seq.random");
