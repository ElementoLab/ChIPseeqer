use lib qw(/home/olly/PERL_MODULES);
use Sets;

die "please specify an output dir ..\n" if (!$ARGV[2]);

#
#  load graph
#
open IN, $ARGV[0];

while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    
    $H{ $a[0] }{ $a[1] } = 1 if (!defined($H{ $a[0] }{ $a[1] }) && !defined($H{ $a[1] }{ $a[0] }));
    

}
close IN;




open IN, $ARGV[1];


while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    my $t = shift @a;

    my $a_ref = Sets::removeDuplicates(\@a);

    @a = @$a_ref;

    my $n = scalar(@a);

    open OUT, ">$ARGV[2]/$t.txt";
    #print "$t\n";

    for (my $i=0; $i<$n-1; $i++) {
	for (my $j=$i+1; $j<$n; $j++) {
	    
	    if (defined($H{ $a[$i] }{ $a[$j] }) || defined($H{ $a[$j] }{ $a[$i] })) {
		print OUT "$a[$i]\t$a[$j]\n";
	    }
	    
	}
    }
    close OUT;
    
}
close IN;

