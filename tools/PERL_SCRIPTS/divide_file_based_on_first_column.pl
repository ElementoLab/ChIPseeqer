use lib qw(/home/olly/PERL_MODULES);
use Sets;


my %H = ();

# read matches
while (my $l = <STDIN>) {
    chomp $l;

    my @a = split /\t/, $l, -1;

    push @{ $H{$a[0]} }, $a[1];
    # get a random position
    
}

foreach my $k (keys(%H)) {
    #open OUT, ">$k.txt";
    my @S = sort { $a <=> $b } @{ $H{$k} };
    foreach my $r (@S) {
	print "$k\t$r\n";
    }
    #close OUT;
}
