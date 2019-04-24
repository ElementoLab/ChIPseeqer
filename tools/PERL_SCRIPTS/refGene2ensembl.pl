while (my $l = <STDIN>) {
    chomp $l;
    
    my @a = split /\t/, $l, -1;

    my $name    = $a[0];
    my $chr     = $a[1];
    my $strand  = $a[2]; $strand = $strand eq "+"?1:-1;
    my $exstart = $a[8];
    my $exend   = $a[9];

    my @a1      = split /\,/, $exstart, -1; pop @a1;
    my @a2      = split /\,/, $exend  , -1; pop @a2;

    my $n1      = scalar(@a1);
    my $n2      = scalar(@a2);

    for (my $i=0; $i<$n1; $i++) {
	print "$chr\t$name\t$a1[$i]\t$a2[$i]\t$strand\n";
    }
}
