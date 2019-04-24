
while (my $l = <STDIN>) {
    chomp $l;
    if ($l =~ /^\ /) {
	$l =~ s/\ \ //;
	$l1 = $l;
    } else {
	my @a = split /\t/, $l;
	print "$l1\t$a[0]\n";
    }

   
}
