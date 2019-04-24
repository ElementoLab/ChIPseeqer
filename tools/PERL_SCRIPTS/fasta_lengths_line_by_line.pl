my $n  = undef;
my $s  = 0;
my $nb = 0;
while ($l = <STDIN>) {
    chomp $l;
    if ($l =~ /^>/) {
	if ($nb > 0) {
	    print "$n\t$s\n";
	}
	$s = 0;
	$n = substr($l, 1);
	$nb ++;
    } else {
	$s += length($l);
    }
}

print "$n\t$s\n";
