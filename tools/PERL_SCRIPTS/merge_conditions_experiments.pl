#!/usr/bin/perl
use lib qw(/home/olly/PERL_MODULES);
use Sets;



#my $a_ref_files = Sets::readSet($ARGV[0]);


my %IDX = ();

# up or down regul
my %IDX_U = ();
my %IDX_D = ();


#foreach my $f (@$a_ref_files) {
    
    # read condition file
    my $cond = undef;
    open IN, $ARGV[0] or die "Cannot ..\n";
    my $start = 0;
    while (my $l = <IN>) {
	chomp $l;
	
	if ($l =~ /Cond (\d+)/) {
	    $start = 0; 
	    #print "Cond=$1\n";
	} else {
	    
	    #next if ($start == 20);
	    my @a = split /\t/, $l;
	    $IDX{ $a[0] } ++;
	    $start ++;
	    
	    

	    if ($a[5] eq "U") {
		$IDX_U{$a[0]}++;
	    }
	    
	    if ($a[5] eq "D") {
		$IDX_D{$a[0]}++;
	    }
	    
	    
	}
    }
    close IN;
#}

my @a_tmp = ();

foreach my $k (keys(%IDX)) {
    my @a = ($k, $IDX{$k});
    push @a_tmp, \@a;
}

my @a_tmp_sorted = sort { $b->[1] <=> $a->[1] } @a_tmp;

# output kmer -> conditions
foreach my $r (@a_tmp_sorted) {

 
    my $cu = (defined($IDX_U{$r->[0]})?$IDX_U{$r->[0]}:0);
    my $cd = (defined($IDX_D{$r->[0]})?$IDX_D{$r->[0]}:0);


    print "$r->[0]\t$r->[1]\t$cu\t$cd\n";


}
