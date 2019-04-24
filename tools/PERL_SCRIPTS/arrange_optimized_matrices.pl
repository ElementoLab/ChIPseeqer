#
# takes a WM as input, and rewrite it in the ScanACE way
#

open IN, $ARGV[0];


my $fs = 0;
my $fm = 0;
my $le = undef;
my $tx = "";
while (my $l = <IN>) {
    
    
    
    if ($l =~ /Motif/) {
	$fm = 1;
    } 

    if ($l =~ /\*/) {
	$fs = 1;
    }

    if ($l =~ /[ATGC]+/) {
	$le = length($l) - 1; 
    }

    $tx .= $l;
    
    
}

close IN;


if ($fm == 0) {
    $tx =  "Motif 1\n" . $tx;
}

if ($fs == 0) {
    my $tm = "*" x $le;
    $tx =  $tx . $tm . "\n";
}

open OUT, ">$ARGV[0]";
print OUT $tx;
close OUT;




