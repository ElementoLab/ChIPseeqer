use lib qw(/home/olly/PERL_MODULES);
use Sets;


while (my $l = <STDIN>) {
    chomp $l;
    
    my @a = split /\t/, $l, -1;

    my $n = shift @a;

    my $std = Sets::stddev(\@a);
    my $avg = Sets::average(\@a);

	#print "$std\t$avg\n";

    my @o = ();
    foreach my $r (@a) {
	$r = sprintf("%3.2f", ($r - $avg)/$std);
        push @o, $r;
        #$_ = sprintf("%3.2f", (defined($_)?$_:0)); #log($_));

    }
    
    print "$n\t"; print join("\t", @o); print "\n";
    
}
