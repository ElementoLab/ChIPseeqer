package Categories;
use lib qw(/home/olly/PERL_MODULES);
use Database;
use Hypergeom;
use Sets;
use Table;
use strict;


sub new {
    my $self  = {};
   

    my %t1 = ();
    my %t2 = (); 
    my %t3 = ();
    $self->{NUMBERS}     = \%t1;
    $self->{ANNOTATIONS} = \%t2;
    $self->{CATEGORIES}  = \%t3;
    $self->{NBCATEGORIES}   = 0;
    $self->{BONFERRONI}     = 1;
    $self->{ORF_SET}     = undef;
    $self->{VERBOSE}     = 1;

    bless($self); 
    return $self;
    
}

sub setResource {
    my ($self, $d) = @_;
    $self->{RESOURCE} = $d;
    die if (! -e "$self->{RESOURCE}/annotations.txt");
    die if (! -e "$self->{RESOURCE}/categories.txt");
    die if (! -e "$self->{RESOURCE}/numbers.txt");

    my $ta = Table->new;
    $ta->loadFile("$self->{RESOURCE}/categories.txt");
    $self->{CATEGORIES} = $ta->getIndexColumnsKV(0,1);

    $ta->loadFile("$self->{RESOURCE}/numbers.txt");
    $self->{NUMBERS} = $ta->getIndexColumnsKV(0,1);
    
    open IN, "$self->{RESOURCE}/annotations.txt";
    while (my $l = <IN>) {
	chomp $l;
	my @a = split /\t/, $l;
	my @b = split /\|/, $a[1];
	$self->{ANNOTATIONS}->{ $a[0] } = \@b; 
    }

    
}


sub setBonferroni {
    my ($self, $i) = @_;
    $self->{BONFERRONI} = $i;
}

#
#  
#
sub setVerbose {
    my ($self, $i) = @_;
    $self->{VERBOSE}    = $i;
}


#
#
#
sub readORFSet {

    my ($self, $s_infile) = @_;     
    $self->{ORF_SET} = Sets::readSet($s_infile);

    if ($self->{VERBOSE} == 1) {
	print "Read " . Sets::size($self->{ORF_SET}) . " ORFs\n";
    }
    
}


#
#    fucking worm genes
#
sub renameWormGenes {
    my ($self) = @_;
    
    foreach my $o (@{$self->{ORF_SET}}) {	
	$o =~ s/\.1$//;	
    }

    Sets::printSet($self->{ORF_SET}) if ($self->{VERBOSE});
}


sub setORFSet {
    my ($self, $a_ref) = @_;
    $self->{ORF_SET} = $a_ref;
    
    Sets::printSet($a_ref) if ($self->{VERBOSE});
}


sub getORFset {
    my ($self) = @_;
    
    return $self->{ORF_SET};
}


sub setTotalNbORFS {
     my ($self, $n) = @_;
     
     $self->{TOTAL_NB_ORFS} = $n;    
}

#
#  compute the union of all annotated genes and the ones provided by the usr
#
sub setTotalNbORFSUnion {
    my ($self, $a_ref) = @_;

    if (scalar(@$a_ref) == 0) {
	die "setTotalNbORFSUnion: please provide a non-empty set ..\n";
    }

    my $cnt = Sets::size(Sets::getUnionSet($a_ref, $self->getAnnotatedORFs()));
    $self->setTotalNbORFS($cnt);
}

# get the categories for a given ORF
sub getCategories {
    
    my ($self, $o) = @_;

    return $self->{ANNOTATIONS}->{$o};
    
}


sub getCategoriesText {
    my ($self, $o) = @_;

    my @a = ();
    foreach my $r (@{$self->{ANNOTATIONS}->{$o}}) {
	push @a, $self->{CATEGORIES}->{$r};
    }
    
    return \@a;
}


sub getAnnotatedORFs {
    my ($self) = @_;

    return keys( %{ $self->{ANNOTATIONS} });
}


#
# 
#
sub setPvalueThreshold {
    
    my ($self, $p) = @_;
     
    $self->{PVALUE_THRESHOLD} = $p;  
    
}


#   
#  get the functional enrichment of an array of ORFS
#
sub getFunctionalContent {
    
    my ($self) = @_;
    
    if (scalar(@{$self->{ORF_SET}}) == 0) {
	die "Please provide a non-empty array\n";
    }
    
    if (!defined($self->{TOTAL_NB_ORFS})) {
	die "Please define the total number of ORFs ..\n";
    }

    my %h_cnt = ();
    my $i_cnt = 0;

    Sets::printSetSep($self->{ORF_SET}, "\t"); <STDIN>;

    #
    #   for each ORF in the set, get the function with which it is annotated 
    #
    foreach my $o (@{$self->{ORF_SET}}) {
	print "OO $o\n";
	Sets::printSetSep($self->{ANNOTATIONS}->{$o}, " ");

	#
	# count the categories associated to a given gene
	#
	foreach my $f (@{ $self->{ANNOTATIONS}->{$o} } ) {
	    $h_cnt{$f} ++;
	}
	
	$i_cnt++;
    }

    Sets::printHash(\%h_cnt);


    #
    # now traverse the non-zero categories and calc a p-value
    #
    my @a_res = ();

    my $s2 = scalar(@{$self->{ORF_SET}});

    reset(%h_cnt);
    while (my ($f,$i) = each(%h_cnt)) {
    	
	#
	# get the number of gene annotated with $f
	#
	my $s1 =  $self->{NUMBERS}->{$f};
	
	
	
	next if ($s1 == 0);
	

	

	
	my $p = Hypergeom::cumhyper($i,$s1,$s2,$self->{TOTAL_NB_ORFS});

	
	my %h_tmp = (
		     P           => $p, 
		     OV          => $i, 
		     S1          => $s1, 
		     S2          => $s2, 
		     CATEGORY          => $self->{CATEGORIES}->{$f}, 
		     ID          => $f, 
		     N           => $self->{TOTAL_NB_ORFS}
		     );
	
	#
	# do nothing if p-vqlue is above threshold
	#
	next if (defined($self->{PVALUE_THRESHOLD}) && ( $p > $self->{PVALUE_THRESHOLD} )); 

	next if (defined($self->{BONFERRONI}) && ($p * $self->{NBCATEGORIES} > 0.05));

	push @a_res, \%h_tmp;
    
    }
    
    my @a_res_bis = sort {$a->{P} <=> $b->{P}} @a_res;

    return \@a_res_bis;
}



1;
