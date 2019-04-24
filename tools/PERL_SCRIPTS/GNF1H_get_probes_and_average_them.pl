#
#  read in gene tto probes mapping
#
use lib qw(/home/olly/PERL_MODULES);
use strict;
use Sets;

my @genes  = ();
my %PROBES = ();
my $cnt = 0;
open IN, $ARGV[0];
while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;

    my $g = shift @a;

    push @genes, $g;

    foreach my $o (@a) {
	push @{ $PROBES { $o } },  $cnt;
    }

    $cnt ++;

}
close IN;


my @EXPRESSION = ();
my @COUNTS     = ();

#
#  go thru expression data, keep only goog probes, and add 
#

open IN, $ARGV[1];
my $l = <IN>; print $l;

my $mmax = undef;
while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;

    my $p = shift @a;

    my $mmax = scalar(@a); 

    if (defined($PROBES{ $p })) {
	
	foreach my $i (@{ $PROBES{ $p }}) {
	    $EXPRESSION[ $i ] = Sets::addArrays($EXPRESSION[ $i ], \@a);
	    $COUNTS    [ $i ] ++;
	}
	
    }

    $cnt ++;
    
}
close IN;


my $n = scalar(@genes);

for (my $i=0; $i<$n; $i++) {
    

    
    foreach my $r (@{ $EXPRESSION[ $i ] }) {
	$r = sprintf("%3.2f", $r / $COUNTS[$i]);
    }

    #my $avg = Sets::average($EXPRESSION[ $i ]);

    #foreach my $r (@{ $EXPRESSION[ $i ] }) {
#	$r = sprintf("%3.2f", Sets::log2( 0.000001 + $r / $avg ));
#    }
    

    print $genes[$i];

    if (scalar(@{ $EXPRESSION[ $i ] }) == 0) {
	print "\t0.0" x 1;
    } else {
	print "\t" . join("\t", @{ $EXPRESSION[ $i ] }); 
    }

    print "\n";
    
}
