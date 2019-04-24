package Mips;
use lib qw(/home/olly/PERL_MODULES);
#use strict;
require 'Sets.pm';

#!/usr/bin/perl

use GDBM_File;
use Hypergeom;
use Getopt::Long;

sub new {
    
    my $self  = {};
    $self->{MIPS_DB_DIR} = "/home/olly/COMPARATIVE_YEAST/ORFS/BLAST_REPORTS/BINDINGMAP/MIPS";

    my %t1 = ();
    my %t2 = (); 
    my %t3 = ();
    
    $self->{DB_NBFUNC_TEXT} = \%t1;
    $self->{DB_NBFUNC_NUMB} = \%t2;
    $self->{DB_ORF_NBFUNC}  = \%t3;

    #print $self->{DB_ORF_NBFUNC};
    
    # takes as input a set of ORFs, and look up the most significant MIPS functions in that 
    tie(%{$self->{DB_NBFUNC_TEXT}}, 'GDBM_File', "$self->{MIPS_DB_DIR}/nbfunc_text.db", &GDBM_READER, 0644);
    tie(%{$self->{DB_NBFUNC_NUMB}}, 'GDBM_File', "$self->{MIPS_DB_DIR}/nbfunc_numb.db", &GDBM_READER, 0644);
    tie(%{$self->{DB_ORF_NBFUNC}},  'GDBM_File', "$self->{MIPS_DB_DIR}/orf_nbfunc.db",  &GDBM_READER, 0644);
    
    #print $self->{DB_NBFUNC_NUMB};
    #print "\n";

    $self->{ORF_SET} = [];
    
    $self->{TOTALNBORFS} = 6200; #5538;  # Kellis et al number
    
    $self->{PVALUE_THRESHOLD} = undef; 

    $self->{BONFERRONI} = scalar(keys(%{$self->{DB_NBFUNC_NUMB}}));
    
    bless($self);
    return $self;
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
	$self->setPvalueThreshold(0.05 /  $self->{BONFERRONI_T});
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
# get an array MIPS categories for one specific gene
#
sub getMIPScategory {
    my ($self, $orf, $type) = @_;

    # get a set of functions
    my $s = $self->{DB_ORF_NBFUNC}->{$orf};

    # get an array of functions
    my @a = split /\|/, $s;

    my @a_tmp = ();

    #print keys(%{$self->{DB_NBFUNC_TEXT}});

    # count the functions
    foreach my $f (@a) {
	#print "$f\n";
	
	next if ($f == 214);
	
	#print "toto=" . $self->{DB_NBFUNC_TEXT}->{$f};
	if (defined($type)) { 
	    push @a_tmp, $f;
	} else {
	    push @a_tmp, $self->{DB_NBFUNC_TEXT}->{$f};
	}
	
    }

    return \@a_tmp;
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

	my $a_ref_tmp  = $self->getMIPScategory($o, 2);
	
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

	
	my %h_tmp = (PVALUE => $p, OVERLAP => $i, TOTAL => $s2, TEXT => $self->{DB_NBFUNC_TEXT}->{$f}, NUM => $f);

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
