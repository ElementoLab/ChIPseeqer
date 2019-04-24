#  INPUT : single gene info, 
#
#
use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;


my $ta = Table->new;

# load gene info
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

# create a list
my @A = ();
my $i = 0;
foreach my $r (@$a_ref) {
    my @a_tmp = ($i, $r->[0]);
    push @A, \@a_tmp;
    $i++;
}

# sort 
@A = sort { $b->[1] <=> $a->[1] } @A;

# put indices in an array
my %H = ();
for (my $i=0; $i<$ARGV[1]; $i++) {
    #print $A[$i]->[1] . "\n";
    $H{$A[$i]->[0]} = 1;
}




# load the gene2gene matrix
open IN, $ARGV[2];

my %MAT = ();

#print join("\t", keys(%H)); print "\n";
my $i = 0;
my %ADJ = ();
while (my $l = <IN>) {
    #print "i=$i\n"; 
    if (!defined($H{$i})) {
	$i++;
	next;
    }
    chomp $l;
    my @a = split /\t/, $l;

    my $n = scalar(@a);

    for (my $j=0; $j<$n; $j++) {
	
	if (defined($H{ $j })) {
	    
	    if ($a[$j] > 0.4) {
		#print $a[$j]; print "\t";
		push @{ $ADJ{ $i } }, $j; 
	    }
	} 
	
    }

    #print "\n";
    

    $i++;
}

close IN;



#
#  load in the gene description
#
$ta->loadFile($ARGV[3]);
my $a_ref_names = $ta->getColumn(1);

my %INDEX = ();
foreach my $r (keys(%ADJ)) {
    
    my $n1 = ($a_ref_names->[$r] ne "inf"?$a_ref_names->[$r]:$r);

    #print "$a_ref_names->[$r]\t$r";
    
    foreach my $v ( @{ $ADJ{$r} } ) {
	
	my $n2 = ($a_ref_names->[$v] ne "inf"?$a_ref_names->[$v]:$v);

	next if (defined($INDEX{$n1}{$n2}) || defined($INDEX{$n2}{$n1}) );

	print "$n1 pp $n2\n";
	
	$INDEX{$n1}{$n2} = 1;
	
    }
    
    #print "\n";
}

my @B = keys(%ADJ);
foreach my $r (@B) {
    $r = ($a_ref_names->[$r] ne "inf"?$a_ref_names->[$r]:$r);

}


Sets::writeSet(\@B, "list_genes_osprey.txt");
