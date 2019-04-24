use Hypergeom;
use Getopt::Long;



$s_orthofile = undef;
GetOptions ('file1=s'     => \$s_file1,
            'file2=s'     => \$s_file2,
	    'orth=s'      => \$s_orthofile,
	    'n=s'      => \$n);


if (!$n) {
    $n = 6000;
}

#print @ARGV;
if (length($ARGV[0]) == 0) {
    #die "Usage : overlap_simple.pl --file1= --file2= --n=\n";
	
}

$h_ref_index = getIndex($s_orthofile) if ($s_orthofile);
$n = scalar(keys(%$h_ref_index)) if ($s_orthofile);
#print keys(%$h_ref_index);


print "$s_file1\t$s_file2\t";
print poverlap($s_file1, $s_file2, $n);
print "\n";




sub poverlap {
    
    my ($s_file1, $s_file2, $n) = @_;

    my %h_there = ();

    my $s1 = 0;

    open OUT, ">overlap.txt";
    
    open IN, $s_file1;
    while (my $s = <IN>) {
	chomp $s;

	next if ($s_orthofile && ($h_ref_index->{$s} != 1));

	$h_there{$s} = 1;

	$s1++;
    }
    close IN;


    open IN, $s_file2;
    $s2 = 0;
    $i  = 0;
    while (my $s = <IN>) {
	chomp $s;

	next if ($s_orthofile && ($h_ref_index->{$s} != 1));	

	$i++ if ($h_there{$s} == 1);

	print OUT "$s\n" if ($h_there{$s} == 1);

	$s2++;
    }
    close IN;
    close OUT;

    print "$i\t$s1\t$s2\t$n\t";
    return Hypergeom::cumhyper($i,$s1,$s2,$n);
    
    
}



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
