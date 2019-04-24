use lib qw(/home/elemento/PERL_MODULES);

use Sets;

my $h_ref = Sets::getIncreasingIndex($ARGV[0]);


open IN, $ARGV[1];
my $l = <IN>;
print "$l";
while (my $l = <IN>) {
    chomp $l;

    my @a = split /\t/, $l, -1;
    
    if (defined($h_ref->{ $a[0] })) {
	#$a[0] = $h_ref->{ $a[0] } . "_" . $a[0];
	print join("\t", @a) . "\n";
    }

}
close IN;


