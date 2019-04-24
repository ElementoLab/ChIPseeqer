#!/usr/bin/perl
use lib qw(/home/olly/PERL_MODULES);
use Sets;

my $k = length($ARGV[0]);

my $a_ref_files = Sets::getFiles("./*.$k" . "mers.txt");

foreach my $r (@$a_ref_files) {
    open IN, $r;
    my $cond = undef;
    my $txt  = undef;
    while (my $l = <IN>) {
	chomp $l;
	
	if ($l =~ /(^Cond.+)\).+over/) {
	    #print $1;
	    $cond = $1;
	    $cnt  = 0;
	} else {
	    my @a = split /\t/, $l;
	    $cnt ++;
	    if (($a[0] eq $ARGV[0]) || ($a[0] eq Sets::getComplement($ARGV[0]))) {
		$txt .= "$cond\t$l\t$cnt\n";
	    }
	}
    }
    
    if ($txt) {
	print "$r\n$txt";;
    }
    

    close IN;
    
}
