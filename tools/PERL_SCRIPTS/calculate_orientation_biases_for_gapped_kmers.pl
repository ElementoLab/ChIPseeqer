use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $f1 = $ARGV[1]; # "../../../DATA/ELE_D_300.seq";
my $f2 = $ARGV[2]; # "../../../DATA/BRI_D_300.seq";
my $ng = $ARGV[3]; # 11292;

foreach my $r (@$a_ref) {
    
    my $cs1 = $r->[3];
    
    my $cre = Sets::getComplement($r->[0]);
    $cre =~ s/N/\./g;
    
    my $ex = `/home/olly/PROGRAMS/FASTCOMPARE/recompare -re \"$cre\" -fasta1 $f1 -fasta2 $f2 -nbgenes $ng -twostrand 0`;

    my ($cs2) = $ex =~ /overlap=(\d+)/;

    print "$r->[0]\t$cs1\t$cre\t$cs2\n";


}

