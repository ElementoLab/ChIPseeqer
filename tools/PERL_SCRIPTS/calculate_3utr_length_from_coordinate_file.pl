while (my $l = <STDIN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;

    if ($a[4] == -1) {
	$l = $a[2] - $a[5];
    } else {
	$l = $a[6] - $a[3];
    }
    
    if ($l < 0) {
	die "problem, len is negative\n";
    }

    print "$a[0]\t$l\n";
}
