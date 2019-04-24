#
#  take as input a set of targets and a microRNA sequence
#    align each target to the miRNA using a modified MATRIX
#    calculate free energy score 
#


use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;
use Fasta;
use Vienna;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $vi = Vienna->new;


my $mir = $ARGV[1];
$mir =~ s/u/t/g;
$mir = Sets::getComplement(uc($mir));

my $linker = 'X' x 8;

while (my $a_ref = $fa->nextSeq()) {

    my ($n, $s) = @$a_ref;

    if ($ARGV[2] == 1) {

	#
	# create a file containing both sequence
	#
	
	open OUT, ">toto1.seq";
	print OUT ">$n\n$s\n";
	close OUT;
	
	open OUT, ">toto2.seq";
	print OUT ">mir\n$mir\n";
	close OUT;
	
	system("needle toto1.seq toto2.seq -datafile ../test.mat -outfile toto.out -gapextend 8 -gapopen 2 > /dev/null 2>&1");
	
	#system("cat toto.out");
	
	my $txt = Sets::file2txt("toto.out");
	
	my ($score) = $txt =~ /\# Score\: ([\d\.]+)/;
	
	print "S=$score\n";
    }

    
    if ($ARGV[2] == 2) {
	
	my $srna = lc($s);
	$srna =~ s/t/u/g;
	#$vi->setRNASeq($ARGV[1] . $linker . $srna);
	
	$vi->setRNASeq($srna . $linker . $ARGV[1]);
	
	
	my $a_ref_folds = $vi->fold;
	
	my $thefold  = shift @$a_ref_folds;
	
	#print "E=$thefold->{ENERGY} ($thefold->{FOLD})\n";
	print "$thefold->{ENERGY}\n"; # ($thefold->{FOLD})\n";

    }

    #<STDIN>;
}



#my $ta = Table->new;
#$ta->loadFile($ARGV[0]);


