use lib qw(/home/olly/PERL_MODULES);

use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref_map = $ta->getIndexKV(0,1);


open IN, $ARGV[1];

my %SAVE = ();

while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    my $t = shift @a;

    
    print "$t";
    foreach my $p (@a) {
	if (defined($h_ref_map->{$p})) {
	    $h_ref_map->{$p} =~ s/\//\t/g;
	    print "\t$h_ref_map->{$p}";
	}
    }
    
    print "\n";
    
}
close IN;
