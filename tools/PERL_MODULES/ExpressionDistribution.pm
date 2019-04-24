package ExpressionDistribution;
use lib qw(/home/olly/PERL_MODULES);

use Sets;
require "GnuPlot.pm";

sub new {
    my $self  = {};
    $self->{SETS} = [];
    $self->{DIR_COR} = "/home/olly/COMPARATIVE_YEAST/ORFS/BLAST_REPORTS/BINDINGMAP/SCRIPTS";
    
    $self->{EXP_DATA} = undef;
    $self->{BKG_COR} = undef;

    $self->{FILENAMES} = [];

    $self->{TITLES} = [];

    bless($self);
    return $self;
}


sub setExpressionDataset {
    my ($self, $s_dataset)  = @_;
    
    $self->{EXP_DATA} = $s_dataset;


}


sub setBackgroundCorrelations {
    my ($self, $s_dataset)  = @_;
    
    $self->{BKG_COR} = $s_dataset;

}

#
# add an array of sets (files)
#
sub addSets {
    my ($self, $a_ref)  = @_;

    push @{ $self->{SETS} }, $a_ref;

}


sub addTitles {
    my ($self, $a_ref)  = @_;
    
    push @{ $self->{TITLES} }, $a_ref;

}




sub addOutputFileName {
    
    my ($self, $file)  = @_;

    push @{ $self->{FILENAMES} }, $file;
    
}

#
# calcule les correlations
#
sub calcCorrelations {

    my ($self)  = @_;

    if (!defined($self->{EXP_DATA})) {
	die "Please define an expression dataset\n";
    }

    if (scalar(@{$self->{FILENAMES}}) == 0) {
	die "Please define at least one file name\n";
    }
    
    # temp file
    my $s_tmpcorr = "/tmp/tmp.corr";
    open TMP, ">$s_tmpcorr";
    foreach $r1 (@{ $self->{SETS} } ) {
	foreach $r2 (@$r1 ) {
	    print TMP "$r2\n";
	}
    }
    close TMP;

    my $s_todo = undef;
    my $s_todo = "$self->{DIR_COR}/correlation -sets $s_tmpcorr -exp $self->{EXP_DATA} -suffix .corr";
    #print "Executing : $s_todo\n" if ($s_verbose eq 'T') ;
    system(($s_todo)) == 0 or die "Cannot execute $s_todo\n"; 

    my $j = 0;
    foreach $r1 (@{ $self->{SETS} } ) {

	# for each set of sets, output the distributions
	my @a_sets = ();
	
	foreach $r2 ( @$r1 ) {
	    
	    my $a_ref = Sets::readSet("$r2.corr");

	    print "Warning, size of the set is 0\n" if (scalar(@$a_ref) == 0);
	    
	    push @a_sets, $a_ref;
	    
	}

	# get the number of bins
	my $i_nbbins = 2 * int(log(scalar(@{$a_sets[0]})));
	 
	my @a_dists = ();
	
	my $gp = GnuPlot->new;
	
	
	my $i = 0;
	my $s_args = "";
	foreach $r2 ( @$r1 ) {
	    
	    my $a_ref = Sets::getDistribution($a_sets[$i], $i_nbbins);

	    Sets::writeDistribution($a_ref, "$r2.corr.dist");
	 
	    system("cat $r2.corr.dist");
   
	    $gp->addFile("$r2.corr.dist");

	    $i++;
	}

	$gp->setTitles($self->{TITLES}->[$j]);
	

	if (defined($self->{BKG_COR})) {

	    my $a_ref1 = Sets::readSet($self->{BKG_COR});
	    my $a_ref2 = Sets::getDistribution($a_ref1, $i_nbbins);
	    
	    Sets::writeDistribution($a_ref2, "$self->{BKG_COR}.dist");
	    
	    $gp->addFile("$self->{BKG_COR}.dist");

	    $gp->addTitle("Background");
	    
	}
	
	
	
	$gp->plot($self->{FILENAMES}->[$j]);;
	
	$j++;
    }
    
    
}




1;
