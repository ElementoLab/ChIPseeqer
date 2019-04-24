package COG;
use lib qw(/home/olly/PERL_MODULES);
#use strict;
require 'Sets.pm';

#!/usr/bin/perl

use GDBM_File;
use Hypergeom;
use Getopt::Long;

sub new {
    
    my $self  = {};
    $self->{COG_DB_FILE}    = "/home/olly/DATA/BACTERIA/COG/gene_to_categories_Ecoli.txt";
    $self->{COG_NAMES_FILE} = "/home/olly/DATA/BACTERIA/COG/cog_names.txt";

    
    
    my %t1 = ();
    my %t2 = (); 
    my %t3 = ();
    $self->{COG_NAMES}      = \%t1;
    #$self->{DB_NBFUNC_TEXT} = \%t1;
    #$self->{DB_NBFUNC_NUMB} = \%t2;
    #$self->{DB_ORF_NBFUNC}  = \%t3;

    #print $self->{DB_ORF_NBFUNC};
    
    # takes as input a set of ORFs, and look up the most significant MIPS functions in that 
    #$self->{DB_NBFUNC_TEXT}; 
    #$self->{DB_NBFUNC_NUMB}}, 'GDBM_File', "$self->{MIPS_DB_DIR}/nbfunc_numb.db", &GDBM_READER, 0644);
    #$self->{DB_ORF_NBFUNC}  
    
    #print $self->{DB_NBFUNC_NUMB};
    #print "\n";

    bless($self);

    $self->readDataFile($self->{COG_DB_FILE});

    
    $self->readCOGNames($self->{COG_NAMES_FILE});

    $self->{ORF_SET}          = [];    
    $self->{TOTALNBORFS}      = undef; #6200; #5538;  # Kellis et al number
    $self->{PVALUE_THRESHOLD} = undef; 
    $self->{BONFERRONI}       = 26*31;
    
    
    return $self;
}



sub setSpecies {
    my ($self, $s) = @_;
    $self->{COG_DB_FILE}    = "/home/olly/DATA/BACTERIA/COG/gene_to_categories_$s.txt";
    $self->readDataFile($self->{COG_DB_FILE});
}



sub readCOGNames {
    my ($self, $f) = @_;
    
    open IN, $f;
    while (my $l = <IN>) {
	chomp $l;
	my @a = split /\t/, $l;
	$self->{COG_NAMES}->{ $a[0] } = $a[1]; 
    }
    close IN;

}


sub readDataFile {
    my ($self, $f) = @_;
    
    open IN, $f;
    while (my $l = <IN>) {
	chomp $l;
	my @a = split /\t/, $l;
	my $g = shift @a;
	$self->{DB_ORF_NBFUNC}->{ $g } = \@a; 
	
	foreach my $c (@a) {
	    $self->{DB_NBFUNC_NUMB}->{ $c } ++;
	}
    }
    close IN;

}

sub readORFset {
    my ($self, $s_infile) = @_;
    $self->{ORF_SET} = readSet($s_infile);
}


#
# 
# 
sub setORFset {
    my ($self, $a_ref) = @_;
    $self->{ORF_SET} = $a_ref;
}


#
#
sub setTotalNbORFS {
     my ($self, $n) = @_;
     
     $self->{TOTALNBORFS} = $n;    
}

sub setBonferroni {
    
    my ($self, $p) = @_;
    
    if ($p == 1) {
	$self->setPvalueThreshold(0.05 /  $self->{BONFERRONI});
    }
}
#
# 
#
sub setPvalueThreshold {
    
    my ($self, $p) = @_;
     
    $self->{PVALUE_THRESHOLD} = $p;  
    
}


#
# get the functional enrichment of an array of ORFS
#
sub getFunctionalContent {
    
    my ($self) = @_;
    
    if (scalar(@{$self->{ORF_SET}}) == 0) {
	die "Please provide a non-empty array\n";
    }
    
    my %h_cnt = ();
    my $i_cnt = 0;
    
    foreach my $o (@{$self->{ORF_SET}}) {

	my $a_ref_tmp  = $self->{DB_ORF_NBFUNC}->{$o};
	
	# count the functions
	foreach my $f (@$a_ref_tmp) {
	    $h_cnt{$f} ++;
	}
	
	$i_cnt++;
    }
    # now traverse the non-zero categories and calc a p-value
    
    my @res = ();

    my $s2 = scalar(@{$self->{ORF_SET}});
    
    reset(%h_cnt);
    while (my ($f,$i) = each(%h_cnt)) {
    
	# get the number of gene annotated with 4f
	my $s1 = $self->{DB_NBFUNC_NUMB}->{$f};

	#print "$i,$s1,$s2,$self->{TOTALNBORFS}\n";

	my $p = Hypergeom::cumhyper($i,$s1,$s2,$self->{TOTALNBORFS});

	
	my %h_tmp = (PVALUE => $p, OVERLAP => $i, TOTAL => $s2, TEXT => $self->{COG_NAMES}->{$f}, NUM => $f);

	# do nothing if p-vqlue is above threshold
	next if (defined($self->{PVALUE_THRESHOLD}) && ( $p > $self->{PVALUE_THRESHOLD} )); 
	
	push @res, \%h_tmp;
    
    }

    
    my @res_bis = sort {$a->{PVALUE} <=> $b->{PVALUE}} @res;

    return \@res_bis;
}


#
# get functional content for a specific CATEGORY
#
sub getSpecificEnrichment {
    
    my ($self, $cat) = @_;

    my $tmp;

    if (defined($self->{PVALUE_THRESHOLD})) {
	$tmp = $self->{PVALUE_THRESHOLD};
	
	$self->{PVALUE_THRESHOLD} = undef;
	
    } 

    foreach my $r (@{$self->getFunctionalContent}) {
	
	if ($r->{NUM} == $cat) {
	    $self->{PVALUE_THRESHOLD} = $tmp;
	    return $r; 
	}
    }
    
    $self->{PVALUE_THRESHOLD} = $tmp;    

    return ();
}


1;
