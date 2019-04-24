use strict;
use Bio::Seq;
use lib qw(/home/olly/PERL_MODULES);
use Sequence;
use Sets;
use DataFiles;


if (scalar(@ARGV) == 0) {
    print "Usage : THISPRG predictions1 predictions2\n";
    exit;
}

my $h_ref1 = my_readEnsemblOrthologs($ARGV[0]);

my $h_ref2 = my_readEnsemblOrthologs($ARGV[1]);

my $a_ref  = my_getReciprocalOrthologs($h_ref1, $h_ref2);


foreach my $r (@$a_ref) {
    print "$r->[0]\t$r->[1]\n";
}


#
#  read an Ensembl ortholog prediction
#    returns an index GENENAME -> ARRAY
#    only keeps the best ortholgs !
sub my_readEnsemblOrthologs {
    
    my ($file) = @_;
    
    my %DATA = ();

    open IN1, $file or die "cannot read ortholog file ..\n";
    while (my $l = <IN1>) {
	chomp $l;
	next if ($l =~ /Chromosome/);

	my @a_tmp = split /\t/, $l;

	# next if ortholog already there and score does not improve
	next if (defined($DATA{$a_tmp[0]} && ($DATA{$a_tmp[0]}->[2] > $a_tmp[2]))); 

	$DATA{$a_tmp[0]} = \@a_tmp;
    }

    close IN1;


    # return 1 gene => best ortholog 
    return \%DATA;

}


#  
#
sub my_getReciprocalOrthologs {
    my ($h_ref1, $h_ref2) = @_;

    my @a = ();
    
    foreach my $k1 (keys(%$h_ref1)) {

	
	# get the best ortholog in species 2
	my $bo2 = $h_ref1->{$k1}->[1];
	
	# get the best ortholog in species 1
	my $bo1 = $h_ref2->{$bo2}->[1];
	
	if ($k1 eq $bo1) {
	    my @a_tmp  = ($k1, $bo2);
	    push @a, \@a_tmp;
	}
	
    }
    
    return \@a;
}
