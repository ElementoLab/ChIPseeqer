# input gene2gene, gene@phenotype, threshold1, threshold2

use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;



my $a_ref = Sets::readSet($ARGV[0]);
my $i = 0;

my @a = @$a_ref;
my @a = sort { $b <=> $a } @a;
my $max = $a[100];

#print "$max\n";

foreach my $r (@$a_ref) {
    if ($r > $max) {
	push @genes, $i;
	$H[ $i ] = 1;
    }
    $i++;
}

my $ta = Table->new;
$ta->loadFile("gene_descriptions.txt");
my $a_ref_names = $ta->getColumn(1);

my %EDGES = ();

open IN, $ARGV[1];
my $i = 0;
while (my $l = <IN>) {
    if (defined($H[$i])) {
	chomp $l;
	my @a = split /\t/, $l;
	my $j = 0;
	my @b = ();
	foreach my $r (@a) {
	    if (defined($H[$j])) {
		push @b, $r;

		#print "$a_ref_names->[$i] $a_ref_names->[$j] $r\n";

		#if ($r > 0.37) {

		 #   my ($g1, $g2) = sort ($i, $j);
		    
		 #   $EDGES{ "$g1" . "_$g2" } = "$a_ref_names->[$g1]_$g1 -- $a_ref_names->[$g2]_$g2;\n"; 
		    
		#}
		
	    }
	    $j++;
	}
	print $a_ref_names->[$i] . "\t"; print join("\t", @b); print "\n";
    }
    $i++;
}
close IN;


#open OUT, ">graph.txt";
#print OUT "graph G {\n";
#foreach my $k (keys(%EDGES)) {
#    print OUT $EDGES{ $k };
#}
#print OUT "}\n";
#close OUT;

