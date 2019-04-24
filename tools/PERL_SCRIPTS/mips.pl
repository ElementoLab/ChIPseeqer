#!/usr/bin/perl

use GDBM_File;
use Hypergeom;
use Getopt::Long;

$s_orthofile = undef;
$s_verbose = "F";
GetOptions (
            
	    'verbose=s'       => \$s_verbose,
	    'set=s'       => \$s_set,
	    'orth=s'      => \$s_orthofile
	    );

$dir = "/home/olly/COMPARATIVE_YEAST/ORFS/BLAST_REPORTS/BINDINGMAP/MIPS";

# takes as input a set of ORFs, and look up the most significant MIPS functions in that 
tie(my(%db_nbfunc_text),'GDBM_File',"$dir/nbfunc_text.db", &GDBM_WRCREAT, 0644);
tie(my(%db_nbfunc_numb),'GDBM_File',"$dir/nbfunc_numb.db", &GDBM_WRCREAT, 0644);
tie(my(%db_orf_nbfunc), 'GDBM_File',"$dir/orf_nbfunc.db", &GDBM_WRCREAT, 0644);



$h_ref_index = getIndex($s_orthofile) if ($s_orthofile);


$s_infile = $s_set;

$a_ref_orfs = readSet($s_infile);

$n = 5538;

# s2
$s2 = scalar(@$a_ref_orfs);

#print "Input set contains $s2 orfs\n"; 

# count the number of ORFs in each MIPS category
%h_cnt = ();
my $i_cnt = 0;
foreach $o (@$a_ref_orfs) {

    #print "$o was assigned function " . $db_orf_nbfunc{$o} . "\n";
    
    next if ($s_orthofile && ($h_ref_index->{$o} != 1));
    

    # get a set of functions
    my $s = $db_orf_nbfunc{$o};
    
    # get an array of functions
    my @a = split /\|/, $s;

    # count the functions
    foreach $f (@a) {
	$h_cnt{$f} ++;
    }

    $i_cnt++;
    
}

$s2 = $i_cnt if $s_orthofile;


# now traverse the non-zero categories and calc a p-value

$i_bonferroni = scalar(keys(%h_cnt));

my @res = ();

reset(%h_cnt);
while (my ($f,$i) = each(%h_cnt)) {
    
    # get the number of 
    my $s1 = $db_nbfunc_numb{$f};

    # next if (($s1 == 1) && ($i == 1));

    my $p = Hypergeom::cumhyper($i,$s1,$s2,$n);

    # $p = $p * $i_bonferroni;

    #print  $p . "\t" . $i . "\t" . $s2 . "\t" . $db_nbfunc_text{$f} . "\n"; # if ($p < 0.001);

    my @a_tmp = ($p, $i, $s2, $db_nbfunc_text{$f}, $f);
    
    push @res, \@a_tmp;
    
}


@res_bis = sort {$a->[0] <=> $b->[0]} @res;


foreach my $r (@res_bis) {
    
    print $r->[0] . "\t" . $r->[1] . "\t" . $r->[2] . "\t" . $r->[3] . "\n";
    
    if ($s_verbose eq "T") {
	foreach my $o (@$a_ref_orfs) {
	    # get a set of functions
	    my $s = $db_orf_nbfunc{$o};
	    
	    # get an array of functions
	    my @a = split /\|/, $s;
	    
	    if (in_array($r->[4], @a)) {
		$g = `yeast_genename.pl $o`; chomp $g;
		print "$o $g\n";
	    }
	}
    }
    
}

sub removeDuplicates {
    my ($a_ref) = @_;
    
    my @b = ();
    
    # inverted inex
    my %ix = ();
    foreach my $aa (@$a_ref) {
	push @b, $aa if (!$ix{$aa});
	$ix{$aa} = 1;
    }

    #print join("\n", @b) . "\n";
    
    return \@b;
    
}



sub in_array() {
    my $val = shift(@_);
 
    foreach $elem(@_) {
        if($val == $elem) {
            return 1;
        }
    }
    return 0;
}


# get the overlap set between two sets
sub getOverlapSet {
    my ($a_ref1, $a_ref2) = @_;
    

    my @a_overlap = ();

    my %h_there = ();

    foreach $k (@$a_ref1) {
	$h_there{$k} = 1;
    }

    foreach $k (@$a_ref2) {
	push @a_overlap, $k if ($h_there{$k} == 1);
    }

    return \@a_overlap;
}



# read a set from disk
sub readSet {
    my ($file) = @_;
    
    open IN, $file;
    my @a_tmp = ();
    while (my $l = <IN>) {
	chomp $l;
	push @a_tmp, $l if ($l ne "");
    }
    close IN;

    my $r = removeDuplicates(\@a_tmp);

    return $r;
}



 


untie(%db_nbfunc_text);
untie(%db_nbfunc_numb);
untie(%db_orf_nbfunc);



sub getIndex {
    my ($f) = @_;

    my %i = ();
    
    open IN, $f;
    while (my $l = <IN>) {
	chomp $l;

	$i{$l} = 1;
    }
    close IN;

    return \%i;
}

